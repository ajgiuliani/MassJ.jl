"""
Chromatographic peak analysis on ion-current traces ([`IonCurrent`](@ref)):
peak detection with boundaries, baseline-subtracted integration, and the
standard chromatography shape/quality metrics (FWHM, asymmetry, USP tailing,
theoretical plates, statistical moments, S/N).

Works on any `IonCurrent` axis (retention time, drift time, compensation
voltage). Best run on a smoothed, baseline-corrected trace — see
[`smooth`](@ref) and [`baseline_correction`](@ref).
"""

export ChromPeak, chrom_peaks, chrom_peak

"""
    struct ChromPeak

A detected chromatographic peak and its metrics. All abscissa quantities are in
the trace's `axis` units (e.g. seconds for `:rt`).

    struct ChromPeak
        axis::Symbol       # :rt, :drift, or :cv
        apex::Float64      # abscissa at the apex
        height::Float64    # baseline-subtracted apex intensity
        left::Float64      # left boundary (abscissa)
        right::Float64     # right boundary (abscissa)
        area::Float64      # baseline-subtracted integrated area
        fwhm::Float64      # full width at half maximum
        asymmetry::Float64 # Aₛ = B/A at 10 % height
        tailing::Float64   # USP tailing factor T = (A+B)/2A at 5 % height
        plates::Float64    # theoretical plates N = 5.54·(apex/FWHM)²
        centroid::Float64  # first moment (intensity-weighted mean abscissa)
        variance::Float64  # second central moment
        skewness::Float64  # standardised third central moment
        snr::Float64       # height / baseline noise
    end
"""
struct ChromPeak
    axis::Symbol
    apex::Float64
    height::Float64
    left::Float64
    right::Float64
    area::Float64
    fwhm::Float64
    asymmetry::Float64
    tailing::Float64
    plates::Float64
    centroid::Float64
    variance::Float64
    skewness::Float64
    snr::Float64
end

function Base.show(io::IO, p::ChromPeak)
    print(io, "ChromPeak(", p.axis, " apex=", round(p.apex, digits = 4),
          ", area=", round(p.area, sigdigits = 4),
          ", height=", round(p.height, sigdigits = 4),
          ", FWHM=", round(p.fwhm, sigdigits = 3),
          ", Aₛ=", round(p.asymmetry, digits = 2),
          ", N=", isnan(p.plates) ? NaN : round(Int, p.plates),
          ", S/N=", round(p.snr, digits = 1), ")")
end

# Robust baseline-noise estimate from successive differences (DER-SNR style):
# for Gaussian noise, σ ≈ MAD(diff)/√2.
function _trace_noise(y::AbstractVector{<:Real})
    length(y) < 3 && return 0.0
    d = abs.(diff(y))
    σ = Statistics.median(d) * 1.4826 / sqrt(2)
    return σ > 0 ? σ : eps()
end

# Analyse one peak over the region arrays `xr`, `yr` (left → right). A local
# linear baseline through the two endpoints is subtracted, then all metrics are
# measured on the baseline-subtracted signal.
function _analyze_peak(xr::AbstractVector{<:Real}, yr::AbstractVector{<:Real},
                       axis::Symbol, noise::Real)
    m = length(xr)
    x1, x2 = xr[1], xr[end]
    y1, y2 = yr[1], yr[end]
    bl = x2 == x1 ? fill(float(y1), m) :
         [y1 + (y2 - y1) * (xr[k] - x1) / (x2 - x1) for k in 1:m]
    yc = collect(float.(yr) .- bl)

    ia     = argmax(yc)
    apex   = xr[ia]
    height = yc[ia]

    halfL = _halfmax_crossing(xr, yc, ia, 0.5height, -1)
    halfR = _halfmax_crossing(xr, yc, ia, 0.5height, +1)
    fwhm  = (isnan(halfL) || isnan(halfR)) ? NaN : halfR - halfL

    a10 = _halfmax_crossing(xr, yc, ia, 0.10height, -1)
    b10 = _halfmax_crossing(xr, yc, ia, 0.10height, +1)
    asym = (isnan(a10) || isnan(b10)) ? NaN : (b10 - apex) / (apex - a10)

    a5 = _halfmax_crossing(xr, yc, ia, 0.05height, -1)
    b5 = _halfmax_crossing(xr, yc, ia, 0.05height, +1)
    tail = (isnan(a5) || isnan(b5)) ? NaN : (b5 - a5) / (2 * (apex - a5))

    area = 0.0
    @inbounds for k in 1:m-1
        area += 0.5 * (yc[k] + yc[k+1]) * (xr[k+1] - xr[k])
    end

    w  = max.(yc, 0.0)
    sw = sum(w)
    if sw > 0
        cen = sum(xr .* w) / sw
        var = sum((xr .- cen) .^ 2 .* w) / sw
        skew = var > 0 ? (sum((xr .- cen) .^ 3 .* w) / sw) / var^1.5 : NaN
    else
        cen, var, skew = apex, NaN, NaN
    end

    plates = (isnan(fwhm) || fwhm <= 0) ? NaN : 5.54 * (apex / fwhm)^2
    snr    = noise > 0 ? height / noise : Inf
    return ChromPeak(axis, apex, height, x1, x2, area, fwhm, asym, tail,
                     plates, cen, var, skew, snr)
end

"""
    chrom_peaks(ic::IonCurrent; snr = 3.0, rel_height = 0.0, noise = nothing) -> Vector{ChromPeak}

Detect peaks in an ion-current trace and return their metrics. Local maxima are
located, boundaries are set valley-to-valley (so adjacent peaks are split at the
valley between them), and each candidate is kept when its baseline-subtracted
`height/noise ≥ snr` and `height ≥ rel_height · max(trace)`. The baseline noise
is estimated robustly from the trace unless given via `noise`.

Run on a smoothed, baseline-corrected trace for best results.
"""
function chrom_peaks(ic::IonCurrent; snr::Real = 3.0, rel_height::Real = 0.0,
                     noise::Union{Real,Nothing} = nothing)
    x = ic.x; y = ic.ic
    n = length(y)
    n < 3 && return ChromPeak[]
    nz   = noise === nothing ? _trace_noise(y) : float(noise)
    ymax = maximum(y)
    peaks = ChromPeak[]
    for i in 2:n-1
        (y[i] > y[i-1] && y[i] >= y[i+1]) || continue          # local maximum
        il = _valley_index(y, i, -1)
        ir = _valley_index(y, i, +1)
        ir - il >= 2 || continue                                # need ≥ 3 points
        p = _analyze_peak(view(x, il:ir), view(y, il:ir), ic.axis, nz)
        (p.snr >= snr && p.height >= rel_height * ymax) || continue
        push!(peaks, p)
    end
    return peaks
end

"""
    chrom_peak(ic::IonCurrent, left::Real, right::Real; noise = nothing) -> ChromPeak

Analyse a single peak within the explicit window `[left, right]` of an
ion-current trace (targeted integration), returning its [`ChromPeak`](@ref)
metrics. The window must contain at least two sample points.
"""
function chrom_peak(ic::IonCurrent, left::Real, right::Real;
                    noise::Union{Real,Nothing} = nothing)
    lo, hi = left <= right ? (left, right) : (right, left)
    il = searchsortedfirst(ic.x, lo)
    ir = searchsortedlast(ic.x, hi)
    il < ir || error("chrom_peak: window [$lo, $hi] contains fewer than 2 points.")
    nz = noise === nothing ? _trace_noise(ic.ic) : float(noise)
    return _analyze_peak(view(ic.x, il:ir), view(ic.ic, il:ir), ic.axis, nz)
end
