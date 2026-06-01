"""
Energy-resolved peak yields from a series of mass spectra.

For each file in a series, the average spectrum is integrated over user-defined m/z
windows. The result is a [`YieldCurve`](@ref) — peak intensities indexed by an
external parameter (photon energy, wavelength, collision energy, …) — that can be
plotted directly or post-processed with [`normalize_tic`](@ref) / [`normalize_flux`](@ref).
"""

# User Interface.
# ---------------

export AbstractPeak, Peak, TargetPeak, YieldCurve,
       yields, integrate_window, normalize_tic, normalize_flux,
       read_peaklist, drop_peaks,
       combine_yields, shift_x, scale_yields, recalibrate_x,
       trim_yields, restrict_x,
       fragment_peaks


const _YIELDS_SUPPORTED_EXT = ("mzml", "mzxml", "mgf", "msp", "imzml", "txt")


# --- Constructors for Peak / TargetPeak (single m/z + tolerance) ------------

"""
    Peak(mz::Real, label::AbstractString; tol = nothing, ppm = nothing) -> Peak
Construct a [`Peak`](@ref) with a fixed window centred on `mz`. Provide exactly
one of `tol` (absolute Δm/z) or `ppm` (parts per million); the resulting window
is `[mz - Δ, mz + Δ]` with `Δ = tol` or `Δ = mz * ppm * 1e-6`.
"""
function Peak(mz::Real, label::AbstractString;
              tol::Union{Real,Nothing} = nothing,
              ppm::Union{Real,Nothing} = nothing)
    Δ = _resolve_tol(mz, tol, ppm, "Peak")
    return Peak(Float64(mz) - Δ, Float64(mz) + Δ, label)
end


"""
    TargetPeak(mz::Real, label::AbstractString;
               tol = nothing, ppm = nothing,
               method::Symbol = :local_max, edges::Real = 0.1) -> TargetPeak
Construct a [`TargetPeak`](@ref) with target `mz`. The search half-width is set
from `tol` (absolute) or `ppm` (parts per million). `method` selects the
per-file resolution algorithm — see [`TargetPeak`](@ref).
"""
function TargetPeak(mz::Real, label::AbstractString;
                    tol::Union{Real,Nothing} = nothing,
                    ppm::Union{Real,Nothing} = nothing,
                    method::Symbol = :local_max,
                    edges::Real    = 0.1)
    Δ = _resolve_tol(mz, tol, ppm, "TargetPeak")
    method ∈ (:local_max, :edges, :centroid) ||
        error("TargetPeak: method must be :local_max, :edges, or :centroid (got :$method)")
    return TargetPeak(Float64(mz), String(label), Δ, method, Float64(edges))
end


function _resolve_tol(mz::Real, tol, ppm, who::String)
    if tol === nothing && ppm === nothing
        error("$who: provide either `tol` (absolute Δm/z) or `ppm`")
    elseif tol !== nothing && ppm !== nothing
        error("$who: provide only one of `tol` or `ppm`")
    end
    return tol === nothing ? Float64(mz) * Float64(ppm) * 1e-6 : Float64(tol)
end


# --- Per-file peak resolution ----------------------------------------------

# Returns (mz1, mz2, found_mz). For a fixed Peak, found_mz = NaN.
_resolve_peak(::MScontainer, ::Any, p::Peak) = (p.mz1, p.mz2, NaN)

function _resolve_peak(spec::MScontainer, centroided, p::TargetPeak)
    if p.method === :centroid
        return _resolve_centroid(centroided, p)
    end
    lo = p.mz - p.tol
    hi = p.mz + p.tol
    idx = findall(x -> lo <= x <= hi, spec.mz)
    if isempty(idx)
        @debug "TargetPeak: no samples in [$lo, $hi] for target $(p.mz) ($(p.label))"
        return (lo, hi, NaN)
    end
    k_rel = argmax(@view spec.int[idx])
    k     = idx[k_rel]
    found = spec.mz[k]
    if p.method === :local_max
        return (found - p.tol, found + p.tol, found)
    else  # :edges
        thresh = p.edges * spec.int[k]
        l = k
        while l > 1 && spec.int[l - 1] > thresh
            l -= 1
        end
        r = k
        while r < length(spec.int) && spec.int[r + 1] > thresh
            r += 1
        end
        return (spec.mz[l], spec.mz[r], found)
    end
end

function _resolve_centroid(centroided, p::TargetPeak)
    centroided === nothing &&
        error("TargetPeak :centroid requires `yields(...; centroid_method=...)`")
    lo = p.mz - p.tol
    hi = p.mz + p.tol
    cmz = centroided.mz
    cit = centroided.int
    idx = findall(x -> lo <= x <= hi, cmz)
    if isempty(idx)
        @debug "TargetPeak :centroid: no centroid in [$lo, $hi] for $(p.mz) ($(p.label))"
        return (lo, hi, NaN)
    end
    k_rel = argmax(@view cit[idx])
    found = cmz[idx[k_rel]]
    return (found - p.tol, found + p.tol, found)
end


_window_of(p::Peak)       = (p.mz1, p.mz2)
_window_of(p::TargetPeak) = (p.mz - p.tol, p.mz + p.tol)

_needs_centroid(peaks) = any(p -> p isa TargetPeak && p.method === :centroid, peaks)


"""
    integrate_window(mz::AbstractVector, int::AbstractVector, mz1::Real, mz2::Real)
    integrate_window(scan::MScontainer, mz1::Real, mz2::Real)
Trapezoidal integration of intensity over the m/z window `[mz1, mz2]`. Returns `0.0`
when fewer than 2 sample points fall inside the window. The order of `mz1`/`mz2` does
not matter; the smaller is taken as the lower bound.
"""
function integrate_window(mz::AbstractVector, intens::AbstractVector, mz1::Real, mz2::Real)
    lo, hi = mz1 <= mz2 ? (mz1, mz2) : (mz2, mz1)
    idx = findall(x -> lo <= x <= hi, mz)
    length(idx) < 2 && return 0.0
    m = @view mz[idx]
    y = @view intens[idx]
    area = 0.0
    @inbounds for k in 1:length(m)-1
        area += 0.5 * (y[k] + y[k+1]) * (m[k+1] - m[k])
    end
    return area
end

integrate_window(scan::MScontainer, mz1::Real, mz2::Real) =
    integrate_window(scan.mz, scan.int, mz1, mz2)


"""
    yields(files::Vector{<:AbstractString}, peaks::Vector{<:AbstractPeak};
           x::AbstractVector{<:Real},
           xlabel::AbstractString = "energy",
           centroid_method::Union{MethodType,Nothing} = nothing)
Build a [`YieldCurve`](@ref) from an explicit list of spectrum files. Each file is
loaded and reduced to a single spectrum with `average(f)`; each
[`AbstractPeak`](@ref) is then resolved against that spectrum (fixed window for
[`Peak`](@ref), located per-file for [`TargetPeak`](@ref)) and integrated.

The external parameter (energy, wavelength, CE…) is supplied either as an
explicit vector `x` (one value per file) **or** as a regular grid
`x = x0 + step*(i-1)` via the `x0` / `step` keywords — matching the convenience
form of [`yields(::AbstractString, ...)`](@ref). Exactly one of the two
conventions must be given.

`centroid_method` is required only when at least one [`TargetPeak`](@ref) uses
`method = :centroid`; it is forwarded to [`MassJ.centroid`](@ref). Defaults to
`MassJ.SNRA(1.0, 100)` in that case.

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

# Examples
```julia-repl
julia> peaks = [Peak(100.5, "A"; tol = 0.5),
                TargetPeak(200.0, "B"; ppm = 5.0, method = :edges)];

julia> yc = yields(["e0.mzML", "e1.mzML"], peaks;
                   x = [3.5, 4.0], xlabel = "photon energy (eV)");

julia> yc = yields(["e0.mzML", "e1.mzML"], peaks;
                   x0 = 3.5, step = 0.5, xlabel = "photon energy (eV)");
```
"""
function yields(files::Vector{<:AbstractString}, peaks::Vector{<:AbstractPeak};
                x::Union{AbstractVector{<:Real},Nothing} = nothing,
                x0::Union{Real,Nothing} = nothing,
                step::Union{Real,Nothing} = nothing,
                xlabel::AbstractString = "energy",
                centroid_method::Union{MethodType,Nothing} = nothing)
    if x !== nothing
        (x0 === nothing && step === nothing) ||
            error("yields: pass either `x` or both `x0` and `step`, not both.")
        length(x) == length(files) ||
            error("yields: length(x) ($(length(x))) != length(files) ($(length(files)))")
        xv = collect(Float64, x)
    elseif x0 !== nothing && step !== nothing
        xv = [Float64(x0) + Float64(step) * (i - 1) for i in 1:length(files)]
    else
        error("yields: provide either `x = [...]` (one value per file) or both `x0` and `step`.")
    end
    nfiles    = length(files)
    npeaks    = length(peaks)
    Y         = Array{Float64}(undef, nfiles, npeaks)
    Y_err     = fill(NaN, nfiles, npeaks)
    found_mz  = fill(NaN, nfiles, npeaks)
    tic       = Vector{Float64}(undef, nfiles)
    tic_err   = fill(NaN, nfiles)

    do_centroid = _needs_centroid(peaks)
    cmethod = do_centroid && centroid_method === nothing ? SNRA(1.0, 100) : centroid_method

    for (i, f) in enumerate(files)
        # Load each file as a vector of scans. The averaged spectrum is used
        # only for *locating* the peak window (stable position); the error on
        # each integrated area is computed from the variance across per-scan
        # integrals, which correctly accounts for inter-bin correlation
        # within a peak.
        scans    = load(f)
        spec     = average(scans)
        cen      = do_centroid ? centroid(spec; method = cmethod) : nothing
        nscans   = length(scans)

        rowtic     = 0.0
        rowvar_acc = 0.0
        any_err    = false
        for (p, peak) in enumerate(peaks)
            lo, hi, found  = _resolve_peak(spec, cen, peak)
            found_mz[i, p] = found

            # Integrate the located window in every scan, then take the
            # standard error of the mean across those per-scan areas.
            areas = Vector{Float64}(undef, nscans)
            @inbounds for k in 1:nscans
                areas[k] = integrate_window(scans[k], lo, hi)
            end
            Y[i, p] = sum(areas) / nscans
            σ = nscans > 1 ? std(areas; corrected = true) / sqrt(nscans) : NaN
            Y_err[i, p] = σ
            rowtic += Y[i, p]
            if isfinite(σ)
                rowvar_acc += σ * σ
                any_err = true
            end
        end
        tic[i]     = rowtic
        tic_err[i] = any_err ? sqrt(rowvar_acc) : NaN
    end
    windows = [_window_of(peak) for peak in peaks]
    labels  = [peak.label        for peak in peaks]
    return YieldCurve(xv, String(xlabel),
                      Y, Y_err, tic, tic_err, found_mz,
                      labels, windows, String.(files), Dict{String,Any}())
end


"""
    yields(dir::AbstractString, peaks::Vector{<:AbstractPeak};
           x0::Real, step::Real,
           xlabel::AbstractString = "energy",
           centroid_method::Union{MethodType,Nothing} = nothing)
Convenience method: list supported spectrum files in `dir` (natural-sort order) and
assign `x = x0 + step*(i-1)` to file `i`.
"""
function yields(dir::AbstractString, peaks::Vector{<:AbstractPeak};
                x0::Real, step::Real,
                xlabel::AbstractString = "energy",
                centroid_method::Union{MethodType,Nothing} = nothing)
    files = list_spectra(dir)
    if isempty(files)
        exts = join(_YIELDS_SUPPORTED_EXT, ", ")
        error("yields: no supported spectra in $dir (extensions: $exts)")
    end
    x = [x0 + step * (i - 1) for i in 1:length(files)]
    return yields(files, peaks; x = x, xlabel = xlabel,
                  centroid_method = centroid_method)
end


# Natural-order sort key: zero-pad digit runs so lexical compare sorts numerically.
_natkey(s::AbstractString) = replace(String(s), r"\d+" => m -> lpad(m, 20, '0'))


"""
    list_spectra(dir::AbstractString) -> Vector{String}
Return paths of all supported spectrum files in `dir`, sorted in natural order
(so "e2.mzML" comes before "e10.mzML").
"""
function list_spectra(dir::AbstractString)
    isdir(dir) || error("list_spectra: directory not found: $dir")
    selected = String[]
    for entry in readdir(dir; join = true)
        isfile(entry) || continue
        ext = lowercase(splitext(entry)[2])
        ext = startswith(ext, ".") ? ext[2:end] : ext
        ext in _YIELDS_SUPPORTED_EXT && push!(selected, entry)
    end
    sort!(selected, by = _natkey)
    return selected
end


"""
    read_peaklist(path::AbstractString;
                  tol::Real = 0.5,
                  ppm::Union{Real,Nothing} = nothing,
                  method::Symbol = :local_max) -> Vector{<:AbstractPeak}
Parse a CSV peak list. The format is auto-detected from the column count:

| Cols | Layout                              | Result                                |
|------|-------------------------------------|---------------------------------------|
| 2    | `mz, label`                         | [`TargetPeak`](@ref) using `tol`/`ppm`/`method` kwargs |
| 3    | `mz1, mz2, label`                   | [`Peak`](@ref) (legacy fixed-window form)             |
| 4    | `mz, tol, method, label`            | [`TargetPeak`](@ref) with per-row `tol` and `method` |

A header row is optional and auto-detected: row 1 is treated as a header when its
first cell is non-numeric. The `tol`, `ppm`, and `method` keywords apply only to
the 2-column form.
"""
function read_peaklist(path::AbstractString;
                       tol::Real = 0.5,
                       ppm::Union{Real,Nothing} = nothing,
                       method::Symbol = :local_max)
    isfile(path) || error("read_peaklist: file not found: $path")
    data  = readdlm(path, ',')
    ncols = size(data, 2)
    ncols in (2, 3, 4) ||
        error("read_peaklist: expected 2, 3, or 4 columns (got $ncols)")
    first_cell = data[1, 1]
    has_header = !(first_cell isa Number) &&
                 tryparse(Float64, strip(string(first_cell))) === nothing
    rows  = has_header ? data[2:end, :] : data
    nrows = size(rows, 1)
    peaks = AbstractPeak[]

    if ncols == 3
        for r in 1:nrows
            push!(peaks, Peak(_tofloat(rows[r, 1]),
                              _tofloat(rows[r, 2]),
                              String(strip(string(rows[r, 3])))))
        end
    elseif ncols == 2
        for r in 1:nrows
            label = String(strip(string(rows[r, 2])))
            mz    = _tofloat(rows[r, 1])
            push!(peaks,
                  ppm === nothing ?
                      TargetPeak(mz, label; tol = tol, method = method) :
                      TargetPeak(mz, label; ppm = ppm, method = method))
        end
    else  # ncols == 4
        for r in 1:nrows
            push!(peaks, TargetPeak(_tofloat(rows[r, 1]),
                                    String(strip(string(rows[r, 4])));
                                    tol    = _tofloat(rows[r, 2]),
                                    method = Symbol(strip(string(rows[r, 3])))))
        end
    end
    return peaks
end

_tofloat(x::Number) = Float64(x)
_tofloat(x)         = parse(Float64, strip(string(x)))

# Lenient: NaN on empty / unparseable / missing.
_try_tofloat(x::Number) = Float64(x)
function _try_tofloat(x)
    s = strip(string(x))
    isempty(s) && return NaN
    p = tryparse(Float64, s)
    return p === nothing ? NaN : p
end


"""
    normalize_tic(yc::YieldCurve) -> YieldCurve
Return a new YieldCurve with each row's peak integrals divided by that row's TIC
(sum of peak integrals across all windows). The `tic` and `tic_err` fields retain
the original raw values. Rows with TIC ≤ 0 are left unchanged.

Errors propagate by the standard division rule
`σ(y/T)² = (1/T)²·σ_y² + (y/T²)²·σ_T²`. The correlation between `y` and `T`
(since `T = Σy`) is ignored — this is the usual first-order approximation.
"""
function normalize_tic(yc::YieldCurve)
    Y       = copy(yc.yields)
    Y_err   = copy(yc.yields_err)
    npeaks  = size(Y, 2)
    for i in 1:size(Y, 1)
        t = yc.tic[i]
        if t > 0
            σ_t = yc.tic_err[i]
            for p in 1:npeaks
                y       = yc.yields[i, p]
                σ_y     = yc.yields_err[i, p]
                Y[i, p] = y / t
                if isfinite(σ_y) && isfinite(σ_t)
                    Y_err[i, p] = sqrt((σ_y / t)^2 + (y * σ_t / (t * t))^2)
                end
            end
        end
    end
    md = copy(yc.metadata)
    md["normalize_tic"] = true
    return YieldCurve(copy(yc.x), yc.xlabel,
                      Y, Y_err, copy(yc.tic), copy(yc.tic_err),
                      copy(yc.found_mz), copy(yc.labels), copy(yc.windows),
                      copy(yc.files), md)
end


"""
    normalize_flux(yc::YieldCurve, flux_file::AbstractString;
                   flux_err_pct::Real = 0.10,
                   skipstart::Int = 0,
                   extrapolate::Symbol = :clamp) -> YieldCurve
Return a new YieldCurve with each row's peak integrals and TIC divided by the photon
flux at that row's `x` value, linearly interpolated from `flux_file`.

`extrapolate` controls behaviour when `yc.x[i]` falls *outside* the flux file's
range:
* `:clamp` (default) — value clamped to the nearest endpoint, a warning is
  emitted.
* `:line` — value linearly extrapolated using the slope of the nearest
  segment (equivalent to `Interpolations.LinearInterpolation(...; extrapolation_bc = Line())`).
  No warning. Same scheme is applied to σ_φ (with `abs(...)` to keep
  uncertainties non-negative).

The flux file has either:

* **2 columns** `x, φ` — the uncertainty on the flux is taken as
  `flux_err_pct * φ` (default 10%); or
* **3 columns** `x, φ, σ_φ` — the third column carries the per-point 1-σ
  uncertainty on the flux and is interpolated alongside `φ`.

Lines starting with `#` are treated as comments and ignored; leading non-numeric
rows are auto-detected and skipped. Use `skipstart = N` to force `N` physical
lines to be discarded from the top of the file before parsing.

Errors propagate by the standard division rule
`σ(y/φ) = (y/φ)·sqrt((σ_y/y)² + (σ_φ/φ)²)`, applied to both the peak yields and
to `tic`. Rows whose interpolated flux is non-positive are left unchanged with a
warning.
"""
function normalize_flux(yc::YieldCurve, flux_file::AbstractString;
                        flux_err_pct::Real  = 0.10,
                        skipstart::Int      = 0,
                        extrapolate::Symbol = :clamp)
    extrapolate ∈ (:clamp, :line) ||
        error("normalize_flux: extrapolate must be :clamp or :line (got :$extrapolate)")
    xf, ff, σf = _read_flux(flux_file; skipstart = skipstart,
                            flux_err_pct = flux_err_pct)
    Y       = copy(yc.yields)
    Y_err   = copy(yc.yields_err)
    tic     = copy(yc.tic)
    tic_err = copy(yc.tic_err)
    npeaks  = size(Y, 2)

    for i in 1:size(Y, 1)
        φ,    in_range  = _interp_linear(xf, ff, yc.x[i]; mode = extrapolate)
        σφ_raw, _       = _interp_linear(xf, σf, yc.x[i]; mode = extrapolate)
        σφ              = abs(σφ_raw)   # extrapolation of σ can go negative
        if !in_range
            @warn "normalize_flux: x=$(yc.x[i]) outside flux range " *
                  "[$(xf[1]), $(xf[end])]; clamped to nearest" flux = φ
        end
        if !(φ > 0)
            @warn "normalize_flux: non-positive flux at x=$(yc.x[i]); " *
                  "skipping division" flux = φ
            continue
        end

        # Peak yields
        for p in 1:npeaks
            y       = yc.yields[i, p]
            σ_y     = yc.yields_err[i, p]
            Y[i, p] = y / φ
            if isfinite(σ_y)
                Y_err[i, p] = sqrt((σ_y / φ)^2 + (y * σφ / (φ * φ))^2)
            elseif y != 0.0
                # propagate flux fraction even when σ_y is unknown — gives at
                # least an estimate of the relative uncertainty
                Y_err[i, p] = abs(y / φ) * abs(σφ / φ)
            end
        end

        # TIC
        t      = yc.tic[i]
        σ_t    = yc.tic_err[i]
        tic[i] = t / φ
        if isfinite(σ_t)
            tic_err[i] = sqrt((σ_t / φ)^2 + (t * σφ / (φ * φ))^2)
        elseif t != 0.0
            tic_err[i] = abs(t / φ) * abs(σφ / φ)
        end
    end
    md = copy(yc.metadata)
    md["normalize_flux"]         = String(flux_file)
    md["normalize_flux_err_pct"] = Float64(flux_err_pct)
    md["normalize_flux_extrap"]  = String(extrapolate)
    return YieldCurve(copy(yc.x), yc.xlabel,
                      Y, Y_err, tic, tic_err,
                      copy(yc.found_mz), copy(yc.labels), copy(yc.windows),
                      copy(yc.files), md)
end


"""
    drop_peaks(yc::YieldCurve, labels) -> YieldCurve
Return a new [`YieldCurve`](@ref) with the peaks whose label is in `labels`
removed. `labels` accepts a single `String` or any iterable of strings; labels
not present in `yc.labels` are silently ignored.

The `tic` field is left unchanged so it still reflects the totals over the
*original* peak set — useful when you want to keep the same TIC reference
after dropping a dominant peak (e.g. for plotting fragment yields without
the precursor swamping the axes).

```julia
plot(drop_peaks(yc, "precursor"))
plot(drop_peaks(yc, ["precursor", "solvent"]))
```
"""
drop_peaks(yc::YieldCurve, label::AbstractString) = drop_peaks(yc, (label,))

function drop_peaks(yc::YieldCurve, labels)
    drop = Set(String.(labels))
    keep = [!(l ∈ drop) for l in yc.labels]
    return YieldCurve(copy(yc.x), yc.xlabel,
                      yc.yields[:, keep], yc.yields_err[:, keep],
                      copy(yc.tic), copy(yc.tic_err),
                      yc.found_mz[:, keep],
                      yc.labels[keep], yc.windows[keep],
                      copy(yc.files), copy(yc.metadata))
end


function _read_flux(path::AbstractString;
                    skipstart::Int   = 0,
                    flux_err_pct::Real = 0.10)
    isfile(path) || error("flux file not found: $path")
    skipstart >= 0 || error("flux file: skipstart must be >= 0 (got $skipstart)")
    raw = readdlm(path; comments = true, comment_char = '#', skipstart = skipstart)
    size(raw, 2) >= 2 || error("flux file must have at least 2 columns (x, flux)")

    nrows = size(raw, 1)
    start = 1
    while start <= nrows
        c = raw[start, 1]
        is_num = (c isa Number) || tryparse(Float64, strip(string(c))) !== nothing
        is_num && break
        start += 1
    end
    start > nrows && error("flux file: no numeric data found in $path")

    rows = raw[start:end, :]
    x = [_tofloat(v) for v in rows[:, 1]]
    y = [_tofloat(v) for v in rows[:, 2]]
    # Column 3 may be entirely absent, partially empty (jagged rows / trailing
    # whitespace), or fully populated. Parse leniently; fall back to
    # `flux_err_pct * |φ|` for rows where the third cell is empty or unparseable.
    pct = Float64(flux_err_pct)
    σ = if size(rows, 2) >= 3
        out = Vector{Float64}(undef, length(y))
        for k in eachindex(y)
            v = _try_tofloat(rows[k, 3])
            out[k] = isfinite(v) ? v : abs(y[k]) * pct
        end
        out
    else
        abs.(y) .* pct
    end
    keep  = .!(isnan.(x) .| isnan.(y) .| isnan.(σ))
    x = x[keep]; y = y[keep]; σ = σ[keep]
    order = sortperm(x)
    return x[order], y[order], σ[order]
end


# Linear interpolation. Boundary handling controlled by `mode`:
#   :clamp (default) — clamps to the nearest endpoint outside [x[1], x[end]]
#   :line            — extends the slope of the nearest segment (Interpolations.jl
#                      `extrapolation_bc = Line()` semantics)
# Returns (value, in_range::Bool). For :line, `in_range` is true even when
# extrapolating, so the caller doesn't emit a clamping warning.
function _interp_linear(x::AbstractVector, y::AbstractVector, xq::Real;
                        mode::Symbol = :clamp)
    n = length(x)
    if xq <= x[1]
        if mode === :line && xq < x[1] && n >= 2
            slope = (y[2] - y[1]) / (x[2] - x[1])
            return y[1] + slope * (xq - x[1]), true
        end
        return y[1], xq == x[1]
    elseif xq >= x[end]
        if mode === :line && xq > x[end] && n >= 2
            slope = (y[end] - y[end - 1]) / (x[end] - x[end - 1])
            return y[end] + slope * (xq - x[end]), true
        end
        return y[end], xq == x[end]
    end
    idx = searchsortedfirst(x, xq)
    x0, x1 = x[idx - 1], x[idx]
    y0, y1 = y[idx - 1], y[idx]
    return y0 + (xq - x0) * (y1 - y0) / (x1 - x0), true
end


"""
    write_csv(yc::YieldCurve, path::AbstractString)
Write a [`YieldCurve`](@ref) to a CSV file: the first row is the header
`<xlabel>, <peak labels…>, TIC`; subsequent rows are one per input file.
"""
function write_csv(yc::YieldCurve, path::AbstractString)
    open(path, "w") do io
        write(io, join(vcat(yc.xlabel, yc.labels, ["TIC"]), ","), "\n")
        for i in 1:length(yc.x)
            row = vcat(yc.x[i], yc.yields[i, :], yc.tic[i])
            write(io, join(row, ","), "\n")
        end
    end
    return path
end


# --- YieldCurve transforms (each returns a new, immutable YieldCurve) --------

# Row-wise TIC and its 1-σ from a (scaled) yields / error matrix — matches the
# convention used when `yields` first builds the curve (TIC = Σ peak areas).
function _recompute_tic(Y::AbstractMatrix, Yerr::AbstractMatrix)
    n = size(Y, 1)
    tic     = vec(sum(Y, dims = 2))
    tic_err = fill(NaN, n)
    @inbounds for i in 1:n
        acc = 0.0; any_err = false
        for j in 1:size(Yerr, 2)
            e = Yerr[i, j]
            if isfinite(e); acc += e * e; any_err = true; end
        end
        tic_err[i] = any_err ? sqrt(acc) : NaN
    end
    return tic, tic_err
end

"""
    combine_yields(ycs::YieldCurve...) -> YieldCurve
Concatenate two or more [`YieldCurve`](@ref)s along the x axis — e.g. stitch
adjacent energy/parameter ranges into one curve. All inputs must share the same
peak columns (identical `labels`); rows are pooled and re-sorted by x.
"""
function combine_yields(ycs::YieldCurve...)
    isempty(ycs) && error("combine_yields: need at least one YieldCurve")
    length(ycs) == 1 && return ycs[1]
    ref = ycs[1]
    for yc in ycs[2:end]
        yc.labels == ref.labels ||
            error("combine_yields: YieldCurves have different peak labels")
    end
    x       = vcat((yc.x          for yc in ycs)...)
    Y       = vcat((yc.yields     for yc in ycs)...)
    Yerr    = vcat((yc.yields_err for yc in ycs)...)
    tic     = vcat((yc.tic        for yc in ycs)...)
    tic_err = vcat((yc.tic_err    for yc in ycs)...)
    fmz     = vcat((yc.found_mz   for yc in ycs)...)
    files   = vcat((yc.files      for yc in ycs)...)
    o = sortperm(x)
    md = Dict{String,Any}("combined_from" => length(ycs))
    return YieldCurve(x[o], ref.xlabel, Y[o, :], Yerr[o, :], tic[o], tic_err[o],
                      fmz[o, :], ref.labels, ref.windows, files[o], md)
end

"""
    shift_x(yc::YieldCurve, Δ::Real) -> YieldCurve
Offset the x axis by `Δ` (everything else unchanged).
"""
function shift_x(yc::YieldCurve, Δ::Real)
    md = copy(yc.metadata)
    md["x_shift"] = get(md, "x_shift", 0.0) + Δ
    return YieldCurve(yc.x .+ Δ, yc.xlabel, yc.yields, yc.yields_err, yc.tic,
                      yc.tic_err, yc.found_mz, yc.labels, yc.windows, yc.files, md)
end

"""
    scale_yields(yc::YieldCurve, factor) -> YieldCurve
Multiply the yields matrix by `factor` — a scalar applied to every peak, or a
length-`npeaks` vector applied per peak. `yields_err` scales by `|factor|`, and
the TIC (and its error) are recomputed from the scaled matrix.
"""
function scale_yields(yc::YieldCurve, factor)
    npeaks = size(yc.yields, 2)
    f = factor isa Real ? fill(Float64(factor), npeaks) :
        (length(factor) == npeaks ? collect(Float64, factor) :
         error("scale_yields: factor must be a scalar or a length-$npeaks vector"))
    Y    = similar(yc.yields)
    Yerr = similar(yc.yields_err)
    for j in 1:npeaks
        @views Y[:, j]    .= yc.yields[:, j]     .* f[j]
        @views Yerr[:, j] .= yc.yields_err[:, j] .* abs(f[j])
    end
    tic, tic_err = _recompute_tic(Y, Yerr)
    md = copy(yc.metadata)
    md["scaled"] = factor isa Real ? Float64(factor) : collect(Float64, factor)
    return YieldCurve(copy(yc.x), yc.xlabel, Y, Yerr, tic, tic_err,
                      copy(yc.found_mz), yc.labels, yc.windows, yc.files, md)
end

"""
    recalibrate_x(yc::YieldCurve, cal::Calibration) -> YieldCurve
    recalibrate_x(yc::YieldCurve, obs_x, ref_x; model=:linear, degree=0) -> YieldCurve
Recalibrate the x axis with a [`Calibration`](@ref) (or fit one on the fly from
`obs_x`/`ref_x`). Rows are re-sorted by the corrected x.
"""
function recalibrate_x(yc::YieldCurve, cal::Calibration)
    newx = cal(yc.x)
    o = sortperm(newx)
    md = copy(yc.metadata)
    md["x_recalibrated"] = string(cal.model)
    return YieldCurve(newx[o], yc.xlabel, yc.yields[o, :], yc.yields_err[o, :],
                      yc.tic[o], yc.tic_err[o], yc.found_mz[o, :], yc.labels,
                      yc.windows, yc.files[o], md)
end

recalibrate_x(yc::YieldCurve, obs_x::AbstractVector{<:Real}, ref_x::AbstractVector{<:Real};
              model::Symbol = :linear, degree::Integer = 0) =
    recalibrate_x(yc, calibrate(obs_x, ref_x; model = model, degree = degree))


# --- Row-wise trimming / restriction of a YieldCurve -------------------------
# (`trim_yields` and `restrict_x` are the row analogue of `drop_peaks`.)

# Build a new YieldCurve keeping the rows whose indices are in `keep` (in order).
function _take_rows(yc::YieldCurve, keep::AbstractVector{<:Integer}, md_note::Pair)
    md = copy(yc.metadata); push!(md, md_note)
    return YieldCurve(yc.x[keep], yc.xlabel,
                      yc.yields[keep, :], yc.yields_err[keep, :],
                      yc.tic[keep], yc.tic_err[keep],
                      yc.found_mz[keep, :], yc.labels, yc.windows,
                      yc.files[keep], md)
end

"""
    trim_yields(yc::YieldCurve; first = 0, last = 0) -> YieldCurve
Drop the first `first` rows and the last `last` rows of a [`YieldCurve`](@ref) —
useful for snipping a startup transient or end-of-scan ringing.

    trim_yields(yc::YieldCurve, drop::AbstractVector{<:Integer}) -> YieldCurve
Drop a list of explicit row indices (1-based into `yc.x`). Out-of-range indices
raise an error; the order of `drop` does not matter.

# Examples
```julia-repl
julia> trim_yields(yc; first = 1, last = 2);    # snip 1 from the start, 2 from the end
julia> trim_yields(yc, [4, 17]);                 # drop two known outliers
```
"""
function trim_yields(yc::YieldCurve; first::Integer = 0, last::Integer = 0)
    first >= 0 && last >= 0 ||
        error("trim_yields: `first` and `last` must be non-negative (got first=$first, last=$last).")
    n = length(yc.x)
    lo = Int(first) + 1
    hi = n - Int(last)
    keep = lo <= hi ? (lo:hi) : (1:0)
    return _take_rows(yc, collect(keep), "trimmed_edges" => (Int(first), Int(last)))
end

function trim_yields(yc::YieldCurve, drop::AbstractVector{<:Integer})
    n = length(yc.x)
    for i in drop
        1 <= i <= n ||
            error("trim_yields: drop index $i is out of range (1:$n).")
    end
    keep = setdiff(1:n, drop)
    return _take_rows(yc, keep, "trimmed_rows" => sort(unique(Int.(drop))))
end

"""
    restrict_x(yc::YieldCurve, xmin::Real, xmax::Real) -> YieldCurve
Keep only the rows whose abscissa `x` lies in `[xmin, xmax]` (inclusive). The
order of `xmin`/`xmax` does not matter; the smaller is used as the lower bound.
"""
function restrict_x(yc::YieldCurve, xmin::Real, xmax::Real)
    lo, hi = xmin <= xmax ? (xmin, xmax) : (xmax, xmin)
    keep = findall(x -> lo <= x <= hi, yc.x)
    return _take_rows(yc, keep, "x_restricted" => (Float64(lo), Float64(hi)))
end


# --- Peptide fragment ions → peak list (bridge to the yields machinery) ------

"""
    fragment_peaks(ions::AbstractVector{FragmentIon}; tol=0.5, ppm=nothing, fixed=true, method=:local_max)
    fragment_peaks(sequence; ions=(:a,:b,:c,:x,:y,:z), charges=1:1, hshifts=(0,), tol=0.5, ppm=nothing, fixed=true, method=:local_max)

Turn peptide fragment ions (see [`fragment_ions`](@ref)) into a peak list for the
[`yields`](@ref) machinery, so spectral/fragmentation yields can be computed.
Each ion becomes a peak centred on its m/z and labelled with the ion label
(e.g. `"b3"`, `"y5(2+)"`, `"w3'"`).

With `fixed = true` (default) each ion becomes a fixed-window [`Peak`](@ref) — a
consistent integration window across a whole spectral series, and quiet when a
fragment is absent. With `fixed = false` each becomes a [`TargetPeak`](@ref) that
is *located* in every spectrum (`method` ∈ `:local_max`, `:edges`, `:centroid`).
Give either `tol` (Da, default 0.5) or `ppm`. The second form builds the ions
from a sequence in one step.

# Examples
```julia-repl
julia> pks = fragment_peaks("PEPTIDE"; ions = (:b, :y, :c, :z), charges = 1:2, ppm = 10);

julia> yc = yields(dir, pks; x0 = 10.0, step = 0.5, xlabel = "collision energy (eV)");
```
"""
function fragment_peaks(ions::AbstractVector{FragmentIon};
                        tol::Real = 0.5, ppm::Union{Real,Nothing} = nothing,
                        fixed::Bool = true, method::Symbol = :local_max)
    peaks = AbstractPeak[]
    for f in ions
        if fixed
            push!(peaks, ppm === nothing ? Peak(f.mz, f.label; tol = tol) :
                                           Peak(f.mz, f.label; ppm = ppm))
        else
            push!(peaks, ppm === nothing ? TargetPeak(f.mz, f.label; tol = tol, method = method) :
                                           TargetPeak(f.mz, f.label; ppm = ppm, method = method))
        end
    end
    return peaks
end

fragment_peaks(sequence::AbstractString;
               ions = (:a, :b, :c, :x, :y, :z), charges = 1:1, hshifts = (0,),
               tol::Real = 0.5, ppm::Union{Real,Nothing} = nothing,
               fixed::Bool = true, method::Symbol = :local_max) =
    fragment_peaks(fragment_ions(sequence; ions = ions, charges = charges, hshifts = hshifts);
                   tol = tol, ppm = ppm, fixed = fixed, method = method)
