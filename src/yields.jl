"""
Energy-resolved peak yields from a series of mass spectra.

For each file in a series, the average spectrum is integrated over user-defined m/z
windows. The result is a [`YieldCurve`](@ref) — peak intensities indexed by an
external parameter (photon energy, wavelength, collision energy, …) — that can be
plotted directly or post-processed with [`normalize_tic`](@ref) / [`normalize_flux`](@ref).
"""

# User Interface.
# ---------------

export Peak, YieldCurve, yields, integrate_window, normalize_tic, normalize_flux, read_peaklist


const _YIELDS_SUPPORTED_EXT = ("mzml", "mzxml", "mgf", "msp", "imzml", "txt")


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
    yields(files::Vector{<:AbstractString}, peaks::Vector{Peak};
           x::AbstractVector{<:Real}, xlabel::AbstractString = "energy")
Build a [`YieldCurve`](@ref) from an explicit list of spectrum files. Each file is
loaded and reduced to a single spectrum with `average(f)`; each [`Peak`](@ref) window
is then integrated. `x` (one value per file) carries the external parameter.

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

# Examples
```julia-repl
julia> peaks = [Peak(100.0, 101.0, "A"), Peak(200.0, 201.0, "B")];

julia> yc = yields(["e0.mzML", "e1.mzML", "e2.mzML"], peaks;
                   x = [3.5, 4.0, 4.5], xlabel = "photon energy (eV)");
```
"""
function yields(files::Vector{<:AbstractString}, peaks::Vector{Peak};
                x::AbstractVector{<:Real},
                xlabel::AbstractString = "energy")
    length(x) == length(files) ||
        error("yields: length(x) ($(length(x))) != length(files) ($(length(files)))")
    nfiles = length(files)
    npeaks = length(peaks)
    Y   = Array{Float64}(undef, nfiles, npeaks)
    tic = Vector{Float64}(undef, nfiles)
    for (i, f) in enumerate(files)
        spec = average(f)
        rowtic = 0.0
        for p in 1:npeaks
            a = integrate_window(spec, peaks[p].mz1, peaks[p].mz2)
            Y[i, p] = a
            rowtic += a
        end
        tic[i] = rowtic
    end
    windows = [(p.mz1, p.mz2) for p in peaks]
    labels  = [p.label       for p in peaks]
    return YieldCurve(collect(Float64, x), String(xlabel), Y, tic,
                      labels, windows, String.(files), Dict{String,Any}())
end


"""
    yields(dir::AbstractString, peaks::Vector{Peak};
           x0::Real, step::Real, xlabel::AbstractString = "energy")
Convenience method: list supported spectrum files in `dir` (natural-sort order) and
assign x = x0 + step*(i-1) to file i.
"""
function yields(dir::AbstractString, peaks::Vector{Peak};
                x0::Real, step::Real,
                xlabel::AbstractString = "energy")
    files = list_spectra(dir)
    if isempty(files)
        exts = join(_YIELDS_SUPPORTED_EXT, ", ")
        error("yields: no supported spectra in $dir (extensions: $exts)")
    end
    x = [x0 + step * (i - 1) for i in 1:length(files)]
    return yields(files, peaks; x = x, xlabel = xlabel)
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
    read_peaklist(path::AbstractString) -> Vector{Peak}
Parse a CSV file with columns `mz1, mz2, label` (with or without a header) into a
`Vector{Peak}`.
"""
function read_peaklist(path::AbstractString)
    isfile(path) || error("read_peaklist: file not found: $path")
    data = readdlm(path, ',')
    size(data, 2) >= 3 ||
        error("read_peaklist: need at least 3 columns (mz1, mz2, label)")
    first = data[1, 1]
    has_header = !(first isa Number) &&
                 tryparse(Float64, strip(string(first))) === nothing
    rows = has_header ? data[2:end, :] : data
    peaks = Peak[]
    for r in 1:size(rows, 1)
        push!(peaks, Peak(_tofloat(rows[r, 1]),
                          _tofloat(rows[r, 2]),
                          String(strip(string(rows[r, 3])))))
    end
    return peaks
end

_tofloat(x::Number) = Float64(x)
_tofloat(x)         = parse(Float64, strip(string(x)))


"""
    normalize_tic(yc::YieldCurve) -> YieldCurve
Return a new YieldCurve with each row's peak integrals divided by that row's TIC
(sum of peak integrals across all windows). The `tic` field retains the original
raw totals. Rows with TIC ≤ 0 are left unchanged.
"""
function normalize_tic(yc::YieldCurve)
    Y = copy(yc.yields)
    for i in 1:size(Y, 1)
        t = yc.tic[i]
        if t > 0
            @views Y[i, :] ./= t
        end
    end
    md = copy(yc.metadata)
    md["normalize_tic"] = true
    return YieldCurve(copy(yc.x), yc.xlabel, Y, copy(yc.tic),
                      copy(yc.labels), copy(yc.windows),
                      copy(yc.files), md)
end


"""
    normalize_flux(yc::YieldCurve, flux_file::AbstractString) -> YieldCurve
Return a new YieldCurve with each row's peak integrals and TIC divided by the photon
flux at that row's `x` value, linearly interpolated from `flux_file`. The flux file
must have 2 header lines followed by 2 whitespace-separated columns (x, flux). Rows
with NaN values are dropped before interpolation; rows whose interpolated flux is
non-positive are left unchanged (a warning is emitted).
"""
function normalize_flux(yc::YieldCurve, flux_file::AbstractString)
    xf, ff = _read_flux(flux_file)
    Y   = copy(yc.yields)
    tic = copy(yc.tic)
    for i in 1:size(Y, 1)
        φ, in_range = _interp_linear(xf, ff, yc.x[i])
        if !in_range
            @warn "normalize_flux: x=$(yc.x[i]) outside flux range " *
                  "[$(xf[1]), $(xf[end])]; clamped to nearest" flux = φ
        end
        if !(φ > 0)
            @warn "normalize_flux: non-positive flux at x=$(yc.x[i]); " *
                  "skipping division" flux = φ
            continue
        end
        @views Y[i, :] ./= φ
        tic[i] /= φ
    end
    md = copy(yc.metadata)
    md["normalize_flux"] = String(flux_file)
    return YieldCurve(copy(yc.x), yc.xlabel, Y, tic,
                      copy(yc.labels), copy(yc.windows),
                      copy(yc.files), md)
end


function _read_flux(path::AbstractString)
    isfile(path) || error("flux file not found: $path")
    raw = readdlm(path; skipstart = 2)
    size(raw, 2) >= 2 || error("flux file must have at least 2 columns (x, flux)")
    x = [_tofloat(v) for v in raw[:, 1]]
    y = [_tofloat(v) for v in raw[:, 2]]
    keep = .!(isnan.(x) .| isnan.(y))
    x = x[keep]; y = y[keep]
    order = sortperm(x)
    return x[order], y[order]
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
