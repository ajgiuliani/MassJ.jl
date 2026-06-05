"""
Processing functions submodule. 
"""

# User Interface.
# ---------------

export smooth, centroid, baseline_correction


"""
    smooth(scan::MScontainer; method::MethodType=SG(5, 9))
Smooth the intensity of the input data and returns a similar structure.
# Examples
```julia-repl
julia> scan = load("test.mzXML")[1];

julia> smoothed_data = smooth(scan)
MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=5.08195e6)
```
"""
function smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
    if method isa MassJ.SG
        return savitzky_golay_filtering(scan, method.order, method.window, method.derivative)
    end  
end

"""
    smooth(scans::AbstractVector{MSscans}; method::MethodType=SG(5, 9, 0))
Smooth the intensity of the input data and returns a similar structure.
# Examples
```julia-repl
julia> scans = load("test.mzXML");

julia> smoothed_data = smooth(scans)
6-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=5.08195e6)
 MSscans(num=2, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=0.7307 min, tic=9727.2  precursor=1255.5)
 MSscans(num=3, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=2.1379 min, tic=11.3032  precursor=902.33)
 MSscans(num=4, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=3.7578 min, tic=4.8084e6)
 MSscans(num=5, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=4.3442 min, tic=12203.5  precursor=1255.5)
 MSscans(num=6, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=5.7689 min, tic=4.84455  precursor=902.33)
```
"""
function smooth(scans::AbstractVector{MSscans}; method::MethodType=SG(5, 9, 0))
    if method isa MassJ.SG
        sm_scans = Vector{MSscans}(undef, 0)
        for el in scans
            push!(sm_scans, savitzky_golay_filtering(el, method.order, method.window, method.derivative))
        end
        return sm_scans
    end  
end



"""
    savitzky_golay_filtering(scan::MassJ.MScontainer, order::Int, window::Int, deriv::Int)
Savitzky-Golay filtering to remove the high frequency noise of int data within the `MSscans` container.
"""
function savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
    y = savitzky_golay(scan.int, order, window, deriv)
    
    basePeakIntensity = ceil(maximum(y))
    basePeakIndex = num2pnt(y, basePeakIntensity)
    basePeakMz = scan.mz[basePeakIndex]
    
    return MSscans(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
end
 

"""
    centroid(scan::MScontainer; method::MethodType=SNRA(1., 100), metrics=false, noise_region=100)
Peak picking on a single spectrum, returning the detected peaks as an
[`MSscans`](@ref). Available algorithms: Signal-to-Noise Ratio Analysis
(`SNRA`), Template Based Peak Detection (`TBPD`), and Continuous Wavelet
Transform (`CWT`).

When `metrics=true`, per-peak descriptors are computed from the input profile and
stored in the returned spectrum's `metadata` under the keys `"peak_fwhm"`,
`"peak_snr"`, `"peak_area"`, and `"peak_resolution"` (one entry per detected
peak, in the same order). `noise_region` is the structuring-element width (in
points) used for the morphological noise floor behind S/N.
"""
function centroid(scan::MScontainer; method::MethodType=SNRA(1., 100),
                  metrics::Bool=false, noise_region::Int=100)
    if method isa TBPD
        ∆mz = 500.0 / method.resolution       # according to mz / ∆mz  = R, we take the value @ m/z 500
        if method.shape == :gauss
            result = tbpd(scan, gauss, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :lorentz
            result = tbpd(scan, lorentz, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :voigt
            result = tbpd(scan, voigt, ∆mz, convert(Float64,method.threshold))
        else
            error("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
        end
    elseif method isa SNRA
        result = snra(scan, method.threshold, method.region)
    elseif method isa CWT
        result = cwt(scan, method.scales, method.threshold, method.min_length)
    else
        error("Unsupported centroiding method $(typeof(method)). Use SNRA, TBPD, or CWT.")
    end

    metrics && _attach_metrics!(result, scan, noise_region)
    return result
end

"""
    centroid(scans::AbstractVector{MSscans}; method::MethodType=SNRA(1., 100) )
Peak picking algorithm taking an array of `MSscans` as input and returning an object of the same type containing the detected peaks. Available algorithm are : Signal to Noise Ratio (SNR) and Template Based Peak Detection (TBPD). Default method is Signal to Noise Ratio Analysis (SNRA), with default threshold = 1.0 and region = 100.
# Examples
```julia-repl
julia> reduced_data = centroid(scans)
6-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS1+, 109 pts m/z=[226.083, 857.333], rt=0.1384 min, tic=521841.4007076025)
 MSscans(num=2, MS2+ CID@18.0eV, 0 pts m/z=[—, —], rt=0.7307 min, tic=0.0  precursor=1255.5)
 MSscans(num=3, MS3+ PQD@35.0eV, 0 pts m/z=[—, —], rt=2.1379 min, tic=0.0  precursor=902.33)
 MSscans(num=4, MS1+, 116 pts m/z=[235.917, 741.083], rt=3.7578 min, tic=468178.80593264103)
 MSscans(num=5, MS2+ CID@18.0eV, 0 pts m/z=[—, —], rt=4.3442 min, tic=0.0  precursor=1255.5)
 MSscans(num=6, MS3+ PQD@35.0eV, 0 pts m/z=[—, —], rt=5.7689 min, tic=0.0  precursor=902.33)
```
"""
function centroid(scans::AbstractVector{MSscans}; method::MethodType=SNRA(1., 100),
                  metrics::Bool=false, noise_region::Int=100)
    return MSscans[centroid(el; method = method, metrics = metrics,
                            noise_region = noise_region) for el in scans]
end

"""
    snra(scan::MScontainer, thres::Real, region::Int)
Signal to Noise Ratio Analysis method returning the m/z and intensity of the peaks detected.
"""
function snra(scan::MScontainer, thres::Real, region::Int)
    # declaration
    peaks_mz = Vector{Float64}(undef,0)
    peaks_int = Vector{Float64}(undef,0)
    peaks_s = Vector{Float64}(undef,0)
    
    noise  = MassJ.opening(scan.int, region);
    SNR = scan.int ./ noise;
        
    maxi = maximum(scan.int)
    # i ranges over interior points only so SNR[i-1] and SNR[i+1] are both
    # in bounds — fixes a BoundsError at the final index when `thres` is
    # small enough that the last point passes the intensity gate.
    for i = 2:length(scan.mz) - 1
        if scan.int[i] >=  maxi * thres / 100.
            if SNR[i] > SNR[i-1] && SNR[i] > SNR[i+1]
                push!(peaks_int, (scan.int[i] - noise[i]))
                push!(peaks_mz, scan.mz[i])
                if !isempty(scan.s)
                    push!(peaks_s, scan.s[i])
                end
            end
        end
    end
    if size(peaks_int, 1) == 0
        basePeakIntensity = NaN
        basePeakMz = NaN
    else
        basePeakIntensity = maximum(peaks_int)
        basePeakMz = peaks_mz[ num2pnt(peaks_int, basePeakIntensity) ]
    end
    
    return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
end

    
"""
    tbpd(scan::MassJ.MScontainer, shape::Symbol,  R::Real, thres::Real)
Template based peak detection algorithm returning the m/z and intensity of the peaks detected.
"""
function tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)   #template based peak detection
    box = num2pnt(scan.mz, scan.mz[1]+0.4) - 1        # taking a box of 0.5 width m/z
    correlation = zeros(length(scan.mz))
    maxi = maximum(scan.int)
    val = 0.0
    for i = 1:1:length(scan.mz)-box
        level = scan.int[i]
        if level >=  maxi * thres / 100. 
            bkg = 0.0
            p0 = [∆mz, scan.mz[i], level, bkg]
            ydata = [model(el, p0) for el in scan.mz[i:i+box]]
            val = Statistics.cor(scan.int[i:i+box], ydata)
        else
            val = 0.0
        end
        if val >= 0.62
            correlation[i] = val
        else
            correlation[i] = 0.0
        end
    end

    peaks_mz = Vector{Float64}(undef,0)
    peaks_int = Vector{Float64}(undef,0)
    peaks_s = Vector{Float64}(undef,0)

    diff_prev = 0.0
    diff      = 0.0
    
    for i = 2:length(correlation)
        if correlation[i] > correlation[i-1] && correlation[i] > correlation[i+1]
            max_value = maximum( scan.int[i:i+2] )
            max_index = num2pnt(scan.int, max_value)
            push!(peaks_mz, scan.mz[max_index])
            push!(peaks_int, scan.int[max_index])
            if !isempty(scan.s)
                push!(peaks_s, scan.s[max_index])
            end
        end
    end
    basePeakIntensity = maximum(peaks_int)
    basePeakMz = peaks_mz[ num2pnt(peaks_int, basePeakIntensity) ]

    return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
end


    
# Ricker ("Mexican hat") wavelet kernel for a given scale `a` (in points),
# sampled on integer offsets out to ±⌈5a⌉.
function _ricker_kernel(a::Real)
    m  = max(1, ceil(Int, 5a))
    js = -m:m
    c  = 2 / (sqrt(3a) * pi^0.25)
    t2 = (js ./ a) .^ 2
    return c .* (1 .- t2) .* exp.(-t2 ./ 2)
end

# One CWT row: cross-correlation of the signal with the Ricker kernel at scale a
# (zero padding at the edges).
function _cwt_row(y::AbstractVector{<:Real}, a::Real)
    ker = _ricker_kernel(a)
    m   = (length(ker) - 1) ÷ 2
    n   = length(y)
    out = zeros(Float64, n)
    @inbounds for i in 1:n
        s = 0.0
        for (kj, j) in enumerate(-m:m)
            idx = i + j
            (1 <= idx <= n) && (s += y[idx] * ker[kj])
        end
        out[i] = s
    end
    return out
end

# Indices at which `v` is a local maximum within a half-window `w` (positive only).
function _local_maxima(v::AbstractVector{<:Real}, w::Int)
    n = length(v); maxima = Int[]
    @inbounds for i in 1:n
        v[i] > 0 || continue
        lo = max(1, i - w); hi = min(n, i + w)
        ismax = true
        for j in lo:hi
            if v[j] > v[i]; ismax = false; break; end
        end
        ismax && push!(maxima, i)
    end
    return maxima
end

# Link local maxima of the CWT matrix into ridge lines, from the coarsest scale
# down to the finest. Each ridge is a vector of (scale_index, position) pairs.
function _cwt_ridges(W::AbstractMatrix{<:Real}, scales)
    n      = size(W, 2)
    order  = sortperm(collect(scales); rev = true)   # coarse → fine
    maxgap = 3
    ridges = Vector{Vector{Tuple{Int,Int}}}()
    open   = Tuple{Int,Int,Int}[]                    # (ridge id, current pos, gap count)
    for si in order
        w    = max(1, round(Int, scales[si]))
        lmax = _local_maxima(view(W, si, :), w)
        used = falses(length(lmax))
        newopen = Tuple{Int,Int,Int}[]
        for (rid, cpos, gaps) in open
            best, bestd = -1, typemax(Int)
            for (li, p) in enumerate(lmax)
                used[li] && continue
                d = abs(p - cpos)
                if d <= w && d < bestd; bestd, best = d, li; end
            end
            if best != -1
                used[best] = true
                push!(ridges[rid], (si, lmax[best]))
                push!(newopen, (rid, lmax[best], 0))
            elseif gaps + 1 <= maxgap
                push!(newopen, (rid, cpos, gaps + 1))
            end
        end
        for (li, p) in enumerate(lmax)
            used[li] && continue
            push!(ridges, [(si, p)])
            push!(newopen, (length(ridges), p, 0))
        end
        open = newopen
    end
    return ridges
end

"""
    cwt(scan::MScontainer, scales, threshold::Real, min_length::Int)
Continuous Wavelet Transform peak detection (Ricker wavelet) after Du, Kibbe &
Lin, *Bioinformatics* **2006**, 22, 2059 (DOI: 10.1093/bioinformatics/btl355).
Builds the CWT of the intensity over `scales`, links local maxima into
ridge lines, and keeps ridges spanning at least `min_length` scales whose
amplitude exceeds `threshold` times the local noise. Returns the detected peaks
as an [`MSscans`](@ref).
"""
function cwt(scan::MScontainer, scales, threshold::Real, min_length::Int)
    y = scan.int
    n = length(y)
    ns = length(scales)
    (n < 3 || ns == 0) && return MSscans(scan.num, scan.rt, 0.0, Float64[], Float64[],
        scan.level, NaN, NaN, scan.precursor, scan.polarity, scan.activationMethod,
        scan.collisionEnergy, Float64[])

    W = Matrix{Float64}(undef, ns, n)
    for (si, a) in enumerate(scales)
        @views W[si, :] .= _cwt_row(y, float(a))
    end

    fine = argmin(collect(scales))          # finest-scale row index
    absfine = abs.(view(W, fine, :))
    minlen  = min_length > 0 ? min_length : max(1, cld(ns, 3))
    halfwin = max(10, n ÷ 50)               # local-noise window half-width

    peaks_mz  = Float64[]; peaks_int = Float64[]; peaks_s = Float64[]
    for r in _cwt_ridges(W, scales)
        length(r) >= minlen || continue
        # finest-scale member gives the best-localised position
        _, j   = findmin([scales[s] for (s, _) in r])
        ipk    = r[j][2]
        strength = maximum(W[s, p] for (s, p) in r)
        lo = max(1, ipk - halfwin); hi = min(n, ipk + halfwin)
        noise = quantile(view(absfine, lo:hi), 0.95)
        snr = noise > 0 ? strength / noise : Inf
        snr >= threshold || continue
        push!(peaks_mz, scan.mz[ipk]); push!(peaks_int, y[ipk])
        isempty(scan.s) || push!(peaks_s, scan.s[ipk])
    end

    if !isempty(peaks_mz)
        o = sortperm(peaks_mz)
        peaks_mz, peaks_int = peaks_mz[o], peaks_int[o]
        isempty(peaks_s) || (peaks_s = peaks_s[o])
        bpi  = maximum(peaks_int)
        bpm  = peaks_mz[num2pnt(peaks_int, bpi)]
    else
        bpi, bpm = NaN, NaN
    end

    return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level,
                   bpm, bpi, scan.precursor, scan.polarity, scan.activationMethod,
                   scan.collisionEnergy, peaks_s)
end


# Walk outward from the apex index `ip` in direction `dir` (±1) until the profile
# drops to `level`, returning the linearly interpolated m/z at the crossing.
# Returns NaN if the profile never reaches `level` before the array edge.
function _halfmax_crossing(mz, int, ip::Int, level::Real, dir::Int)
    n = length(int)
    i = ip
    while true
        j = i + dir
        (j < 1 || j > n) && return NaN
        if int[j] <= level
            denom = int[i] - int[j]
            t = denom == 0 ? 0.0 : (int[i] - level) / denom
            return mz[i] + t * (mz[j] - mz[i])
        end
        i = j
    end
end

# Walk outward from the apex index `ip` to the nearest local minimum (valley) in
# direction `dir` (±1); returns the valley index (or the array edge).
function _valley_index(int, ip::Int, dir::Int)
    n = length(int)
    i = ip
    while true
        j = i + dir
        (j < 1 || j > n) && return i
        int[j] > int[i] && return i
        i = j
    end
end

# Per-peak metrics from the profile spectrum (`mz`, `int`) at the detected peak
# positions `peaks_mz`. Returns (FWHM, S/N, area, resolution) vectors, one entry
# per peak. Half-maximum is taken above the local morphological noise floor;
# area is the valley-to-valley trapezoidal integral of the profile.
function _peak_metrics(mz::AbstractVector{<:Real}, int::AbstractVector{<:Real},
                       peaks_mz::AbstractVector{<:Real}; region::Int = 100)
    np = length(peaks_mz)
    fwhm = fill(NaN, np); snr = fill(NaN, np); area = fill(NaN, np); reso = fill(NaN, np)
    (np == 0 || length(mz) < 2) && return (fwhm, snr, area, reso)

    r = max(1, min(region, length(int) - 1))
    noise = opening(collect(float.(int)), r)
    for (k, pmz) in enumerate(peaks_mz)
        ip = num2pnt(mz, pmz)
        apex = float(int[ip])
        nz = float(noise[ip])
        snr[k] = nz > 0 ? apex / nz : Inf
        half = nz + 0.5 * (apex - nz)
        lo = _halfmax_crossing(mz, int, ip, half, -1)
        hi = _halfmax_crossing(mz, int, ip, half, +1)
        if !isnan(lo) && !isnan(hi) && hi > lo
            fwhm[k] = hi - lo
            reso[k] = pmz / fwhm[k]
        end
        li = _valley_index(int, ip, -1)
        ri = _valley_index(int, ip, +1)
        area[k] = integrate_window(mz, int, mz[li], mz[ri])
    end
    return (fwhm, snr, area, reso)
end

# Compute per-peak metrics from `profile` at the peaks in `result` and store them
# in `result.metadata` (mutating the metadata dict in place).
function _attach_metrics!(result::MSscans, profile::MScontainer, region::Int)
    fwhm, snr, area, reso = _peak_metrics(profile.mz, profile.int, result.mz; region = region)
    result.metadata["peak_fwhm"]       = fwhm
    result.metadata["peak_snr"]        = snr
    result.metadata["peak_area"]       = area
    result.metadata["peak_resolution"] = reso
    return result
end


"""
    baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
Baseline correction takes an `MSscans` object as input and returns an object of the same type with the mass spectra corrected for their base line. Available algorithms are Top Hat, Locally Weighted Error Sum of Squares regression (LOESS) and Iterative Polynomial Smoothing Algorithm (IPSA). The default method is TopHat with a structuring element width of 100 points.
# Examples
```julia-repl
julia> scan = load("test.mzXML")[1];

julia> reduced_data = baseline_correction(scan)
MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=4.596695088878757e6)

julia> baseline_correction(scan, method = MassJ.LOESS(1));      # locally-weighted regression

julia> baseline_correction(scan, method = MassJ.IPSA(51, 100)); # iterative polynomial smoothing
```
"""
function baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
    if method isa TopHat
        return tophat_filter(scan, method.region)
    elseif method isa LOESS
        return loess(scan, method.iter)
    elseif method isa IPSA
        return ipsa(scan, method.width, method.maxiter)
    end
end


"""
    baseline_correction(scans::AbstractVector{MSscans}; method::MethodType=TopHat(100) )
Baseline correction takes a Vector of `MSscans` as input and returns a Vector of `MSscans` with the mass spectra corrected for their base line. Available algorithms are Top Hat, Locally Weighted Error Sum of Squares regression (LOESS) and Iterative Polynomial Smoothing Algorithm (IPSA). The default method is TopHat with a structuring element width of 100 points.
# Examples
```julia-repl
julia> reduced_data = baseline_correction(scans)
6-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=4.596695088878757e6)
 MSscans(num=2, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=0.7307 min, tic=9727.209462129576  precursor=1255.5)
 MSscans(num=3, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=2.1379 min, tic=11.303162056790438  precursor=902.33)
 MSscans(num=4, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=3.7578 min, tic=4.346634123935582e6)
 MSscans(num=5, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=4.3442 min, tic=12203.47988298906  precursor=1255.5)
 MSscans(num=6, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=5.7689 min, tic=4.844552205043209  precursor=902.33)

julia> baseline_correction(scans, method = MassJ.LOESS(1));     # locally-weighted regression

julia> baseline_correction(scans, method = MassJ.IPSA(51, 100)); # iterative polynomial smoothing
```
"""
function baseline_correction(scans::AbstractVector{MSscans}; method::MethodType=TopHat(100) )
    bl_scans = Vector{MSscans}(undef,0)
    for el in scans
        if method isa TopHat
            push!(bl_scans, tophat_filter(el, method.region))            
        elseif method isa LOESS
            push!(bl_scans, loess(el, method.iter))
        elseif method isa IPSA
            push!(bl_scans, ipsa(el, method.width, method.maxiter))
        end       
    end
    return bl_scans
end




"""
    loess(scan::MScontainer, iter::Int )
Method  taking an `MSscans` object as input and returning an object of the same type with the mass spectra without their base line, using the LOESS (Locally Weighted Error Sum of Squares regression).
"""
# Numeric LOESS-baseline kernel on plain (x, y) vectors — shared by the MSscans
# and IonCurrent `baseline_correction` methods. Returns the baseline-subtracted
# signal `res = y - baseline`.
function _loess(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, iter::Int)
    n = length(x)
    r = Int(ceil(n / 2))
    h = [sort(abs.(x .- x[i]))[r] for i = 1:n]
    w = clamp.(abs.((x .- transpose(x)) ./ h), 0.0, 1.0)
    w = (1 .- w.^3).^3
    baseline = zeros(n)
    res   = zeros(n)
    delta = ones(n)
    for _ = 1:iter
        for i = 1:n
            weight = delta .* w[:, i]
            b = [sum(weight .* y), sum(weight .* (y .* x))]
            A = reshape([sum(weight), sum(weight .* x),
                         sum(weight .* x), sum(weight .* x.^2)], 2, 2)
            beta = LinearAlgebra.pinv(A) * b
            baseline[i] = beta[1] + beta[2] * x[i]
        end
        res = y - baseline
        s = Statistics.median(abs.(res))
        delta = clamp.(res ./ (6.0 .* s), -1, 1)
        delta = (1 .- delta.^2).^2
    end
    return res
end

function loess(scan::MScontainer, iter::Int )
    res = _loess(scan.mz, scan.int, iter)
    basePeakIntensity = maximum(res)
    basePeakMz = scan.mz[num2pnt(scan.int, basePeakIntensity)]
    return MSscans(scan.num, scan.rt, sum(res), scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
end


"""
    ipsa(scan::MScontainer, width::Real, maxiter::Int)
Method  taking an `MSscans` object as input and returning an object of the same type with the mass spectra without their base line, using the iterative polynomial smoothing algorithm (IPSA) baseline correction.
"""
# Numeric IPSA-baseline kernel on a plain signal vector `y` — shared by the
# MSscans and IonCurrent `baseline_correction` methods. Returns `res = y - bkg`.
function _ipsa(y::AbstractVector{<:Real}, width::Integer, maxiter::Integer)
    iseven(width) && (width -= 1)
    eps = 1e-07
    input   = zeros(length(y))
    bkg     = savitzky_golay(y, 0, width, 0)
    bkg_old = zeros(length(y))
    res     = y - bkg
    eratio_old = 0.0
    counter = 1
    while true
        for i = 1:length(input)
            input[i] = y[i] < bkg[i] ? y[i] : bkg[i]
        end
        bkg = savitzky_golay(input, 0, width, 0)
        res = y - bkg
        eratio = norm(bkg - bkg_old) / norm(bkg)
        if abs(eratio - eratio_old) < eps
            break
        elseif counter > maxiter
            break
        end
        counter += 1
        eratio_old = eratio
        bkg_old = bkg
    end
    return res
end

function ipsa(scan::MScontainer, width::Real, maxiter::Int)
    res = _ipsa(scan.int, Int(width), maxiter)
    basePeakIntensity = ceil(maximum(res))
    basePeakMz = scan.mz[num2pnt(res, basePeakIntensity)]
    return MSscans(scan.num, scan.rt, scan.tic, scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
end




"""
    tophat_filter(scan::MScontainer, region::Int)
Method taking a MScontainer object as input and returning an object of the same type with the mass spectra without their base line, using the TopHat filtering algorithm.
"""
function tophat_filter(scan::MScontainer, region::Int )
    TIC = sum( tophat(scan.int, region) )
    basePeakIntensity = maximum(tophat(scan.int, region))
    basePeakMz = scan.mz[num2pnt(scan.int,basePeakIntensity)]
    return MSscans(scan.num, scan.rt, TIC, scan.mz, tophat(scan.int, region), scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
end


# --- Ion-current trace (chromatogram / mobilogram / ionogram) processing -----
# These reuse the same numeric kernels as the MSscans methods, applied to the
# trace's intensity `ic` (and abscissa `x` for LOESS), returning a new IonCurrent
# on the same axis.

"""
    smooth(ic::IonCurrent; method::MethodType = SG(5, 9, 0))
Savitzky-Golay smoothing of an ion-current trace (chromatogram, mobilogram, or
ionogram), returning a new [`IonCurrent`](@ref) on the same axis.
"""
function smooth(ic::IonCurrent; method::MethodType = SG(5, 9, 0))
    method isa SG || error("smooth(::IonCurrent): only SG smoothing is supported.")
    y = savitzky_golay(ic.ic, method.order, method.window, method.derivative)
    return IonCurrent(ic.x, y; axis = ic.axis, mobilityType = ic.mobilityType)
end

"""
    baseline_correction(ic::IonCurrent; method::MethodType = TopHat(100))
Baseline / background removal on an ion-current trace, returning a new
[`IonCurrent`](@ref) on the same axis. Methods: [`MassJ.TopHat`](@ref) (default),
[`MassJ.LOESS`](@ref), [`MassJ.IPSA`](@ref).
"""
function baseline_correction(ic::IonCurrent; method::MethodType = TopHat(100))
    y = method isa TopHat ? tophat(ic.ic, method.region) :
        method isa LOESS  ? _loess(ic.x, ic.ic, method.iter) :
        method isa IPSA   ? _ipsa(ic.ic, method.width, method.maxiter) :
        error("baseline_correction(::IonCurrent): use TopHat, LOESS, or IPSA.")
    return IonCurrent(ic.x, y; axis = ic.axis, mobilityType = ic.mobilityType)
end


