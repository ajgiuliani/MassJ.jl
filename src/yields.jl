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
       read_peaklist, drop_peaks


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
        @warn "TargetPeak: no samples in [$lo, $hi] for target $(p.mz) ($(p.label))"
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
        @warn "TargetPeak :centroid: no centroid in [$lo, $hi] for $(p.mz) ($(p.label))"
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


# Internal: trapezoidal integration that also returns the 1-σ uncertainty on the
# area, derived from the per-m/z standard error of the mean (SEM) carried by an
# averaged MSscans. Falls back to (area, NaN) when no error info is available
# (single scan, MSscan, or only one point in the window).
function _integrate_window_with_err(spec::MScontainer, mz1::Real, mz2::Real)
    lo, hi = mz1 <= mz2 ? (mz1, mz2) : (mz2, mz1)
    idx = findall(x -> lo <= x <= hi, spec.mz)
    n = length(idx)
    n < 2 && return (0.0, NaN)

    area = 0.0
    @inbounds for k in 1:n - 1
        i, j  = idx[k], idx[k + 1]
        area += 0.5 * (spec.int[i] + spec.int[j]) * (spec.mz[j] - spec.mz[i])
    end

    # Per-point variance is available only on MSscans with > 1 averaged scan
    if spec isa MSscans && length(spec.num) > 1 && !isempty(spec.s)
        N = length(spec.num)
        var_acc = 0.0
        @inbounds for k in 1:n
            i = idx[k]
            w = if k == 1
                0.5 * (spec.mz[idx[2]] - spec.mz[idx[1]])
            elseif k == n
                0.5 * (spec.mz[idx[n]] - spec.mz[idx[n - 1]])
            else
                0.5 * (spec.mz[idx[k + 1]] - spec.mz[idx[k - 1]])
            end
            sem_i_sq = spec.s[i] / (N * (N - 1))
            var_acc += w * w * sem_i_sq
        end
        return (area, sqrt(max(var_acc, 0.0)))
    end
    return (area, NaN)
end


"""
    yields(files::Vector{<:AbstractString}, peaks::Vector{<:AbstractPeak};
           x::AbstractVector{<:Real},
           xlabel::AbstractString = "energy",
           centroid_method::Union{MethodType,Nothing} = nothing)
Build a [`YieldCurve`](@ref) from an explicit list of spectrum files. Each file is
loaded and reduced to a single spectrum with `average(f)`; each
[`AbstractPeak`](@ref) is then resolved against that spectrum (fixed window for
[`Peak`](@ref), located per-file for [`TargetPeak`](@ref)) and integrated.

`x` (one value per file) carries the external parameter (energy, wavelength, CE…).
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
```
"""
function yields(files::Vector{<:AbstractString}, peaks::Vector{<:AbstractPeak};
                x::AbstractVector{<:Real},
                xlabel::AbstractString = "energy",
                centroid_method::Union{MethodType,Nothing} = nothing)
    length(x) == length(files) ||
        error("yields: length(x) ($(length(x))) != length(files) ($(length(files)))")
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
        spec = average(f)
        cen  = do_centroid ? centroid(spec; method = cmethod) : nothing
        rowtic     = 0.0
        rowvar_acc = 0.0
        any_err    = false
        for (p, peak) in enumerate(peaks)
            lo, hi, found  = _resolve_peak(spec, cen, peak)
            a, σ           = _integrate_window_with_err(spec, lo, hi)
            Y[i, p]        = a
            Y_err[i, p]    = σ
            found_mz[i, p] = found
            rowtic        += a
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
    return YieldCurve(collect(Float64, x), String(xlabel),
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
                   skipstart::Int = 0) -> YieldCurve
Return a new YieldCurve with each row's peak integrals and TIC divided by the photon
flux at that row's `x` value, linearly interpolated from `flux_file`.

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
                        flux_err_pct::Real = 0.10,
                        skipstart::Int     = 0)
    xf, ff, σf = _read_flux(flux_file; skipstart = skipstart,
                            flux_err_pct = flux_err_pct)
    Y       = copy(yc.yields)
    Y_err   = copy(yc.yields_err)
    tic     = copy(yc.tic)
    tic_err = copy(yc.tic_err)
    npeaks  = size(Y, 2)

    for i in 1:size(Y, 1)
        φ,    in_range  = _interp_linear(xf, ff, yc.x[i])
        σφ,   _         = _interp_linear(xf, σf, yc.x[i])
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


# Linear interpolation with nearest-neighbour clamping outside [x[1], x[end]].
# Returns (value, in_range::Bool).
function _interp_linear(x::AbstractVector, y::AbstractVector, xq::Real)
    if xq <= x[1]
        return y[1], xq == x[1]
    elseif xq >= x[end]
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
