"""
Single-run 2-D (m/z × retention-time) feature detection for LC-MS data.

A *feature* is an ion (a narrow m/z) that elutes as a chromatographic peak over a
range of retention times. Detection follows the standard two-step recipe
(à la centWave): build *regions of interest* — mass traces that persist across
consecutive scans at roughly constant m/z — then run chromatographic peak
detection on each trace, reusing [`chrom_peaks`](@ref) so every shape/quality
metric (area, FWHM, S/N, …) comes for free.

Multi-sample alignment / correspondence is deliberately out of scope; this is a
single-run detector.
"""

export Feature, detect_features, feature_table

"""
    struct Feature

A detected 2-D LC-MS feature: a chromatographic peak localised at a specific m/z.

    struct Feature
        mz::Float64        # intensity-weighted m/z of the mass trace
        rt::Float64        # apex retention time
        rt_left::Float64   # left boundary (retention time)
        rt_right::Float64  # right boundary (retention time)
        area::Float64      # baseline-subtracted integrated area
        height::Float64    # baseline-subtracted apex intensity
        fwhm::Float64      # full width at half maximum (retention-time units)
        npoints::Int       # scans contributing within the peak window
        snr::Float64       # height / baseline noise
        mz_min::Float64    # lowest m/z seen in the trace
        mz_max::Float64    # highest m/z seen in the trace
    end
"""
struct Feature
    mz::Float64
    rt::Float64
    rt_left::Float64
    rt_right::Float64
    area::Float64
    height::Float64
    fwhm::Float64
    npoints::Int
    snr::Float64
    mz_min::Float64
    mz_max::Float64
end

function Base.show(io::IO, f::Feature)
    print(io, "Feature(m/z=", round(f.mz, digits = 4),
          ", rt=", round(f.rt, digits = 3),
          " [", round(f.rt_left, digits = 3), ", ", round(f.rt_right, digits = 3), "]",
          ", area=", round(f.area, sigdigits = 4),
          ", height=", round(f.height, sigdigits = 4),
          ", n=", f.npoints,
          ", S/N=", round(f.snr, digits = 1), ")")
end


# --- region of interest (mass trace) being built scan by scan ---------------
mutable struct _ROI
    mzmean::Float64                    # running intensity-weighted mean m/z (match key)
    wsum::Float64                      # Σ intensity (the weight)
    mzmin::Float64
    mzmax::Float64
    first_scan::Int
    last_scan::Int
    pts::Dict{Int,Float64}             # scan index => intensity at this trace
end

_newroi(si::Int, mz::Real, int::Real) =
    _ROI(float(mz), float(int), float(mz), float(mz), si, si,
         Dict{Int,Float64}(si => float(int)))

function _extend!(roi::_ROI, si::Int, mz::Real, int::Real)
    roi.mzmean = (roi.mzmean * roi.wsum + mz * int) / (roi.wsum + int)
    roi.wsum  += int
    roi.mzmin  = min(roi.mzmin, mz)
    roi.mzmax  = max(roi.mzmax, mz)
    roi.last_scan = si
    roi.pts[si] = int
    return roi
end


"""
    detect_features(scans::AbstractVector{MSscans}; ppm = 10.0, mz_tol = nothing,
                    min_scans = 5, max_gap = 1, snr = 3.0, level = 1,
                    centroid_method = nothing) -> Vector{Feature}
    detect_features(filename::AbstractString; kwargs...) -> Vector{Feature}

Detect 2-D (m/z × retention-time) features in a single LC-MS run.

The run is reduced to its MS `level` scans (default `1`), ordered by retention
time, and — when `centroid_method` is given — each scan is centroided first
(otherwise the scans' own peak lists are used). Mass traces are then built: each
peak extends the nearest open trace whose mean m/z is within tolerance (`ppm`, or
an absolute `mz_tol`), starting a new trace otherwise; a trace closes after more
than `max_gap` consecutive scans without a matching peak. Traces spanning at
least `min_scans` scans are passed to [`chrom_peaks`](@ref) (with the given
`snr`) on a zero-filled retention-time trace, and each chromatographic peak found
becomes one [`Feature`](@ref).

Returns the features sorted by `(mz, rt)`. See [`feature_table`](@ref) for a
Tables.jl view.

# Example
```julia
feats = detect_features("lcms_run.mzML"; ppm = 8, min_scans = 6, snr = 5)
```
"""
function detect_features(scans::AbstractVector{MSscans};
                         ppm::Real = 10.0,
                         mz_tol::Union{Real,Nothing} = nothing,
                         min_scans::Integer = 5,
                         max_gap::Integer = 1,
                         snr::Real = 3.0,
                         level::Integer = 1,
                         centroid_method::Union{MethodType,Nothing} = nothing)
    isempty(scans) && return Feature[]

    # 1. select the requested MS level, order by retention time
    ms = [s for s in scans if isempty(s.level) || first(s.level) == level]
    isempty(ms) && (ms = collect(MSscans, scans))     # fallback: no level match
    ms = ms[sortperm([first(s.rt) for s in ms])]
    centroid_method !== nothing &&
        (ms = MSscans[centroid(s; method = centroid_method) for s in ms])
    rts = [first(s.rt) for s in ms]
    nsc = length(ms)

    _tol(mz) = mz_tol === nothing ? mz * float(ppm) * 1e-6 : float(mz_tol)

    # 2. build regions of interest (mass traces)
    open     = _ROI[]
    finished = _ROI[]
    for si in 1:nsc
        mz, intt = ms[si].mz, ms[si].int
        used = falses(length(mz))
        # extend each open trace with the nearest unused peak within tolerance
        for roi in open
            tol = _tol(roi.mzmean)
            best = 0
            bestd = tol
            @inbounds for k in eachindex(mz)
                (used[k] || intt[k] <= 0) && continue
                d = abs(mz[k] - roi.mzmean)
                if d <= bestd
                    bestd = d
                    best = k
                end
            end
            if best != 0
                _extend!(roi, si, mz[best], intt[best])
                used[best] = true
            end
        end
        # close traces idle for more than max_gap scans
        if !isempty(open)
            keep = _ROI[]
            for roi in open
                (si - roi.last_scan > max_gap) ? push!(finished, roi) : push!(keep, roi)
            end
            open = keep
        end
        # start new traces from the still-unused peaks
        @inbounds for k in eachindex(mz)
            (used[k] || intt[k] <= 0) && continue
            push!(open, _newroi(si, mz[k], intt[k]))
        end
    end
    append!(finished, open)

    # 3. chromatographic peak detection per trace → features
    feats = Feature[]
    for roi in finished
        length(roi.pts) >= min_scans || continue
        a, b = roi.first_scan, roi.last_scan
        b - a + 1 >= 3 || continue
        xr = rts[a:b]
        yr = [get(roi.pts, si, 0.0) for si in a:b]
        ic = IonCurrent(collect(Float64, xr), yr; axis = :rt)
        for cp in chrom_peaks(ic; snr = snr)
            npts = count(si -> haskey(roi.pts, si) &&
                               cp.left <= rts[si] <= cp.right, a:b)
            push!(feats, Feature(roi.mzmean, cp.apex, cp.left, cp.right,
                                 cp.area, cp.height, cp.fwhm, npts, cp.snr,
                                 roi.mzmin, roi.mzmax))
        end
    end
    sort!(feats, by = f -> (f.mz, f.rt))
    return feats
end

detect_features(filename::AbstractString; kwargs...) =
    detect_features(load(filename); kwargs...)


"""
    feature_table(features::AbstractVector{Feature}) -> NamedTuple

A [Tables.jl](https://github.com/JuliaData/Tables.jl)-ready column table (one row
per feature) with columns `mz`, `rt`, `rt_left`, `rt_right`, `area`, `height`,
`fwhm`, `npoints`, `snr`, `mz_min`, `mz_max`. Being a `NamedTuple` of vectors it
drops straight into `DataFrame`, `CSV.write`, and the rest of the data ecosystem.

```julia
using DataFrames
df = DataFrame(feature_table(detect_features("run.mzML")))
```
"""
function feature_table(features::AbstractVector{Feature})
    return (mz       = [f.mz       for f in features],
            rt       = [f.rt       for f in features],
            rt_left  = [f.rt_left  for f in features],
            rt_right = [f.rt_right for f in features],
            area     = [f.area     for f in features],
            height   = [f.height   for f in features],
            fwhm     = [f.fwhm     for f in features],
            npoints  = [f.npoints  for f in features],
            snr      = [f.snr      for f in features],
            mz_min   = [f.mz_min   for f in features],
            mz_max   = [f.mz_max   for f in features])
end
