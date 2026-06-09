"""
Energy-resolved peak yields from a series of mass spectra.

For each file in a series, the average spectrum is integrated over user-defined m/z
windows. The result is a [`YieldCurve`](@ref) — peak intensities indexed by an
external parameter (photon energy, wavelength, collision energy, …) — that can be
plotted directly or post-processed with [`normalize_tic`](@ref) / [`normalize_external`](@ref).
"""

# User Interface.
# ---------------

export AbstractPeak, Peak, TargetPeak, YieldCurve,
       yields, integrate_window, normalize_tic, normalize_external, normalize_flux,
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
    TargetPeak(mzs::AbstractVector{<:Real}, label::AbstractString; tol/ppm, …) -> TargetPeak
    TargetPeak(formula::AbstractString, label::AbstractString;
               charge = 1, adduct = "", p_target = 0.99, tol/ppm, …) -> TargetPeak

Construct a [`TargetPeak`](@ref). The search half-width is set from `tol`
(absolute) or `ppm` (parts per million).

- **Single target** `TargetPeak(mz, label; …)` — located per spectrum; `method`
  selects the resolution algorithm (see [`TargetPeak`](@ref)).
- **Explicit cluster** `TargetPeak([mz1, mz2, …], label; …)` — several m/z under
  one label (e.g. isotopologues). Their windows (`mzᵢ ± tol`) are integrated and
  summed into one yield. `ppm`, when given, sets the half-width from the first
  (monoisotopic) m/z.
- **From a formula** `TargetPeak("Nd(NO3)4", label; charge, adduct, p_target, …)`
  — the isotopologue m/z are generated with
  [`isotopic_distribution`](@ref) (keeping isotopologues up to cumulative
  probability `p_target`) and bundled as a cluster. With `adduct = ""` (default)
  the m/z are the neutral isotopologue masses divided by `abs(charge)` (matching
  `isotopic_distribution`); with a non-empty `adduct` (e.g. `"[M+H]+"`) the
  charged m/z are computed with [`adduct_mz`](@ref) and `charge` is taken from the
  adduct.

For the cluster and formula forms `method` selects how the pattern is placed:
`:fixed` (default) integrates fixed windows at the theoretical m/z, while
`:anchor` locates the **anchor** isotopologue in each spectrum and shifts the
whole pattern by the calibration offset (preserving the isotope spacing), falling
back to fixed windows when the anchor is absent. The single-target location
methods (`:local_max` / `:edges` / `:centroid`) are accepted on a cluster for API
symmetry but treated as `:fixed`.

The anchor is chosen at construction. For a **formula** the `anchor` keyword is
`:max` (default — the most-abundant isotopologue, the robust choice for heavy ions
whose monoisotopic peak may be negligible), `:mono` (lowest m/z), or an explicit
m/z. For an **explicit cluster** `anchor` is the monoisotopic m/z by default, or an
explicit m/z that snaps to the nearest cluster member.
"""
function TargetPeak(mz::Real, label::AbstractString;
                    tol::Union{Real,Nothing} = nothing,
                    ppm::Union{Real,Nothing} = nothing,
                    method::Symbol = :local_max,
                    edges::Real    = 0.1)
    Δ = _resolve_tol(mz, tol, ppm, "TargetPeak")
    method ∈ (:local_max, :edges, :centroid) ||
        error("TargetPeak: method must be :local_max, :edges, or :centroid (got :$method)")
    return TargetPeak(Float64(mz), [Float64(mz)], String(label), Δ, method, Float64(edges))
end

function TargetPeak(mzs::AbstractVector{<:Real}, label::AbstractString;
                    tol::Union{Real,Nothing} = nothing,
                    ppm::Union{Real,Nothing} = nothing,
                    method::Symbol = :fixed,
                    edges::Real    = 0.1,
                    anchor::Union{Real,Nothing} = nothing)
    isempty(mzs) && error("TargetPeak: m/z cluster must be non-empty")
    m = sort(Float64.(mzs))
    # The anchor is the m/z that `:anchor` resolution locates in each spectrum.
    # Default = monoisotopic (smallest); an explicit `anchor` snaps to the nearest
    # cluster member (so it is always a real isotopologue position).
    a = anchor === nothing ? m[1] : m[argmin(abs.(m .- Float64(anchor)))]
    Δ = _resolve_tol(a, tol, ppm, "TargetPeak")
    # Clusters resolve with :fixed or :anchor; the single-target location methods
    # are accepted for API symmetry but treated as :fixed (a cluster is not located
    # point-by-point — see `_resolve_windows`).
    method ∈ (:fixed, :anchor, :local_max, :edges, :centroid) ||
        error("TargetPeak: method must be :fixed, :anchor, :local_max, :edges, or :centroid (got :$method)")
    return TargetPeak(a, m, String(label), Δ, method, Float64(edges))
end

function TargetPeak(formula::AbstractString, label::AbstractString;
                    charge::Int = 1, adduct::AbstractString = "",
                    p_target::Real = 0.99,
                    tol::Union{Real,Nothing} = nothing,
                    ppm::Union{Real,Nothing} = nothing,
                    method::Symbol = :fixed,
                    edges::Real    = 0.1,
                    anchor::Union{Symbol,Real} = :max)
    f = MassJ.formula(String(formula))            # dict form avoids the println in the String method
    if isempty(adduct)
        charge == 0 && error("TargetPeak: charge must be non-zero for a bare formula (adduct = \"\")")
        dist = isotopic_distribution(f, p_target; charge = charge)
        mzs  = Float64.(@view dist[2:end, 1])     # m/z column = neutral mass / abs(charge)
    else
        dist    = isotopic_distribution(f, p_target; charge = 1)  # charge = 1 → neutral isotopologue masses
        neutral = Float64.(@view dist[2:end, 1])
        mzs     = [adduct_mz(m, adduct) for m in neutral]         # charge taken from the adduct
    end
    probs = Float64.(@view dist[2:end, 2])
    # Anchor selection. `:max` (default) = most-abundant isotopologue — reliably
    # observable even for heavy ions whose monoisotopic peak is negligible;
    # `:mono` = lowest m/z; a Real picks the nearest isotopologue to that m/z.
    anchor_mz = anchor === :max  ? mzs[argmax(probs)] :
                anchor === :mono ? minimum(mzs)       :
                anchor isa Real  ? Float64(anchor)    :
                error("TargetPeak: anchor must be :max, :mono, or an m/z value (got $anchor)")
    return TargetPeak(mzs, String(label);
                      tol = tol, ppm = ppm, method = method, edges = edges, anchor = anchor_mz)
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

# Returns (windows::Vector{Tuple{Float64,Float64}}, found_mz). One peak may yield
# several integration windows (an isotope cluster); they are summed per scan by
# the caller. `found_mz` is the representative (located or anchor) m/z, or NaN.

_resolve_windows(::MScontainer, ::Any, p::Peak) = ([(p.mz1, p.mz2)], NaN)

function _resolve_windows(spec::MScontainer, centroided, p::TargetPeak)
    if length(p.mzs) > 1
        # Multi-target cluster (e.g. isotopologues). `:anchor` locates the
        # monoisotopic peak and shifts the whole pattern; anything else uses
        # fixed windows at the theoretical m/z.
        return p.method === :anchor ? _resolve_anchored(spec, p) : _resolve_fixed(p)
    end
    lo, hi, found = _resolve_single(spec, centroided, p)
    return ([(lo, hi)], found)
end

# Fixed cluster windows at the theoretical m/z (merged so overlaps count once).
_resolve_fixed(p::TargetPeak) =
    (_merge_windows([(m - p.tol, m + p.tol) for m in p.mzs]), p.mz)

# Anchored cluster: locate the monoisotopic anchor (`p.mz`) in this spectrum,
# then shift every isotopologue window by the calibration offset δ = found − mz
# (so the *theoretical spacing* is preserved while the pattern tracks drift).
# Falls back to the fixed windows, with a debug note, when the anchor is absent.
function _resolve_anchored(spec::MScontainer, p::TargetPeak)
    lo = p.mz - p.tol
    hi = p.mz + p.tol
    idx = findall(x -> lo <= x <= hi, spec.mz)
    if isempty(idx)
        @debug "TargetPeak :anchor: anchor $(p.mz) not found in [$lo, $hi]; using fixed windows ($(p.label))"
        return _resolve_fixed(p)
    end
    found = spec.mz[idx[argmax(@view spec.int[idx])]]
    δ = found - p.mz
    wins = _merge_windows([(m + δ - p.tol, m + δ + p.tol) for m in p.mzs])
    return (wins, found)
end

# Single-target location (the original per-spectrum behaviour). Returns (lo,hi,found).
function _resolve_single(spec::MScontainer, centroided, p::TargetPeak)
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

# Merge overlapping/touching (lo,hi) windows so a shared region is counted once.
function _merge_windows(wins::AbstractVector)
    isempty(wins) && return Tuple{Float64,Float64}[]
    s = sort(wins; by = first)
    out = Tuple{Float64,Float64}[s[1]]
    for (lo, hi) in s[2:end]
        plo, phi = out[end]
        if lo <= phi
            out[end] = (plo, max(phi, hi))
        else
            push!(out, (lo, hi))
        end
    end
    return out
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
_window_of(p::TargetPeak) = length(p.mzs) > 1 ?
    (minimum(p.mzs) - p.tol, maximum(p.mzs) + p.tol) :
    (p.mz - p.tol, p.mz + p.tol)

# Only single-target TargetPeaks locate via centroid; clusters use fixed windows.
_needs_centroid(peaks) =
    any(p -> p isa TargetPeak && length(p.mzs) == 1 && p.method === :centroid, peaks)


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
    yields(files::AbstractVector{<:AbstractString}, peaks::Vector{<:AbstractPeak},
           filters::FilterType...;
           x | x0,step, xlabel = "energy",
           centroid_method = nothing)
Build a [`YieldCurve`](@ref) from an explicit list of spectrum files. Each file is
loaded, **filtered by any [`FilterType`](@ref) arguments** (e.g. `Level(1)`,
`RT(10,20)`, `Activation_Method("CID")`), reduced to a single spectrum with
`average`, and each [`AbstractPeak`](@ref) is resolved against that spectrum
(fixed window for [`Peak`](@ref), located per-file for [`TargetPeak`](@ref)) and
integrated. The integrated-area error is the SEM across the per-scan integrals of
the *filtered* scans.

The external parameter (energy, wavelength, CE…) is supplied either as an
explicit vector `x` (one value per file) **or** as a regular grid
`x = x0 + step*(i-1)` via the `x0` / `step` keywords. Exactly one convention.

`centroid_method` is required only when at least one [`TargetPeak`](@ref) uses
`method = :centroid`; it is forwarded to [`MassJ.centroid`](@ref) (default
`MassJ.SNRA(1.0, 100)`). For full control over per-point preprocessing, pass a
pre-built series instead — see [`yields(::AbstractVector{<:AbstractVector{MSscans}}, ...)`](@ref).

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

# Examples
```julia-repl
julia> peaks = [Peak(100.5, "A"; tol = 0.5),
                TargetPeak(200.0, "B"; ppm = 5.0, method = :edges)];

julia> yc = yields(["e0.mzML", "e1.mzML"], peaks;
                   x = [3.5, 4.0], xlabel = "photon energy (eV)");

julia> yc = yields(["e0.mzML", "e1.mzML"], peaks, Level(1), RT(5.0, 15.0);
                   x0 = 3.5, step = 0.5);          # only MS1 scans in 5–15 min
```
"""
# Resolve the external-parameter vector from either `x` or `x0`/`step` for `n` points.
function _resolve_x(x, x0, step, n::Integer)
    if x !== nothing
        (x0 === nothing && step === nothing) ||
            error("yields: pass either `x` or both `x0` and `step`, not both.")
        length(x) == n ||
            error("yields: length(x) ($(length(x))) != number of points ($n)")
        return collect(Float64, x)
    elseif x0 !== nothing && step !== nothing
        return [Float64(x0) + Float64(step) * (i - 1) for i in 1:n]
    else
        error("yields: provide either `x = [...]` (one value per point) or both `x0` and `step`.")
    end
end

# Apply composed FilterType predicates to a scan list (identity when no filters).
function _filter_scans(scans::AbstractVector{MSscans}, filters)
    isempty(filters) && return scans
    pred = compose_predicates(scans, filters)
    return [s for s in scans if pred(s)]
end

# Core: integrate the peaks over a *series of scan-lists* (one list per point).
# The averaged spectrum locates each window; the area error is the SEM across the
# per-scan integrals, which accounts for inter-bin correlation within a peak.
function _yields_core(series::AbstractVector{<:AbstractVector{MSscans}},
                      peaks::Vector{<:AbstractPeak}, xv::Vector{Float64},
                      xlabel::AbstractString, point_labels::Vector{String},
                      centroid_method)
    npts     = length(series)
    npeaks   = length(peaks)
    Y        = zeros(Float64, npts, npeaks)
    Y_err    = fill(NaN, npts, npeaks)
    found_mz = fill(NaN, npts, npeaks)
    tic      = zeros(Float64, npts)
    tic_err  = fill(NaN, npts)

    do_centroid = _needs_centroid(peaks)
    cmethod = do_centroid && centroid_method === nothing ? SNRA(1.0, 100) : centroid_method

    for (i, scans) in enumerate(series)
        if isempty(scans)
            @debug "yields: no scans for point $i ($(point_labels[i])); row left at zero."
            continue
        end
        spec   = average(scans)
        cen    = do_centroid ? centroid(spec; method = cmethod) : nothing
        nscans = length(scans)

        rowtic     = 0.0
        rowvar_acc = 0.0
        any_err    = false
        for (p, peak) in enumerate(peaks)
            wins, found    = _resolve_windows(spec, cen, peak)
            found_mz[i, p] = found

            # Integrate every (merged) sub-window in each scan and SUM them within
            # the scan, then take the SEM across the per-scan totals — keeping the
            # pattern-level error correct when the windows co-vary (isotopes).
            areas = Vector{Float64}(undef, nscans)
            @inbounds for k in 1:nscans
                a = 0.0
                for (lo, hi) in wins
                    a += integrate_window(scans[k], lo, hi)
                end
                areas[k] = a
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
                      labels, windows, point_labels, Dict{String,Any}())
end

function yields(files::AbstractVector{<:AbstractString}, peaks::Vector{<:AbstractPeak},
                filters::FilterType...;
                x::Union{AbstractVector{<:Real},Nothing} = nothing,
                x0::Union{Real,Nothing} = nothing,
                step::Union{Real,Nothing} = nothing,
                xlabel::AbstractString = "energy",
                centroid_method::Union{MethodType,Nothing} = nothing)
    xv     = _resolve_x(x, x0, step, length(files))
    series = [_filter_scans(load(f), filters) for f in files]
    return _yields_core(series, peaks, xv, xlabel, String.(files), centroid_method)
end

"""
    yields(series::AbstractVector{<:AbstractVector{MSscans}}, peaks, filters::FilterType...;
           x | x0,step, xlabel, labels, centroid_method) -> YieldCurve

Build a [`YieldCurve`](@ref) from a pre-built *series of scan-lists* — one
`Vector{MSscans}` per external-parameter point. This is the composable form: the
caller does any preprocessing (filtering, smoothing, baseline correction,
calibration…) before handing the scans over, and the per-scan SEM is preserved.
Optional [`FilterType`](@ref) arguments are applied to each list as well. Name the
points with `labels` (defaults to `point_1`, `point_2`, …).

```julia
series = [load(f) |> s -> smooth.(s) for f in files]   # any per-point prep
yc = yields(series, peaks, Level(1); x = energies)
```
"""
function yields(series::AbstractVector{<:AbstractVector{MSscans}}, peaks::Vector{<:AbstractPeak},
                filters::FilterType...;
                x::Union{AbstractVector{<:Real},Nothing} = nothing,
                x0::Union{Real,Nothing} = nothing,
                step::Union{Real,Nothing} = nothing,
                xlabel::AbstractString = "energy",
                labels::Union{AbstractVector{<:AbstractString},Nothing} = nothing,
                centroid_method::Union{MethodType,Nothing} = nothing)
    xv  = _resolve_x(x, x0, step, length(series))
    pts = [_filter_scans(collect(MSscans, s), filters) for s in series]
    plabels = labels === nothing ? ["point_$i" for i in 1:length(series)] : String.(labels)
    length(plabels) == length(series) ||
        error("yields: length(labels) ($(length(plabels))) != number of points ($(length(series)))")
    return _yields_core(pts, peaks, xv, xlabel, plabels, centroid_method)
end

"""
    yields(dir::AbstractString, peaks::Vector{<:AbstractPeak}, filters::FilterType...;
           x0, step, x, xlabel, centroid_method)
Convenience method: list supported spectrum files in `dir` (natural-sort order) and
assign `x = x0 + step*(i-1)` to file `i` (or pass an explicit `x`). [`FilterType`](@ref)
arguments are applied to each file's scans before averaging/integration.
"""
function yields(dir::AbstractString, peaks::Vector{<:AbstractPeak}, filters::FilterType...;
                x0::Union{Real,Nothing} = nothing,
                step::Union{Real,Nothing} = nothing,
                x::Union{AbstractVector{<:Real},Nothing} = nothing,
                xlabel::AbstractString = "energy",
                type = nothing,
                centroid_method::Union{MethodType,Nothing} = nothing)
    files = list_spectra(dir; type = type)
    if isempty(files)
        exts = join(_YIELDS_SUPPORTED_EXT, ", ")
        error("yields: no supported spectra in $dir (extensions: $exts)")
    end
    return yields(files, peaks, filters...;
                  x = x, x0 = x0, step = step, xlabel = xlabel,
                  centroid_method = centroid_method)
end


# Natural-order sort key: zero-pad digit runs so lexical compare sorts numerically.
_natkey(s::AbstractString) = replace(String(s), r"\d+" => m -> lpad(m, 20, '0'))


"""
    list_spectra(dir::AbstractString) -> Vector{String}
Return paths of all supported spectrum files in `dir`, sorted in natural order
(so "e2.mzML" comes before "e10.mzML").
"""
function list_spectra(dir::AbstractString; type = nothing)
    isdir(dir) || error("list_spectra: directory not found: $dir")
    allowed = _allowed_exts(type)
    selected = String[]
    for entry in readdir(dir; join = true)
        isfile(entry) || continue
        _ext_of(entry) in allowed && push!(selected, entry)
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
    normalize_external(yc::YieldCurve, external_file::AbstractString;
                       err_pct::Real = 0.10,
                       skipstart::Int = 0,
                       extrapolate::Symbol = :clamp) -> YieldCurve
Return a new YieldCurve with each row's peak integrals and TIC divided by an
**external quantity** (one not contained in the mass spectra) interpolated to that
row's `x` value from `external_file`. The quantity can be anything the yields
should be normalised against: photon flux, laser power or pulse energy, ion-source
current, detector efficiency, transmission, acquisition time — any per-`x`
reference with its own (optional) uncertainty.

`extrapolate` controls behaviour when `yc.x[i]` falls *outside* the file's range:
* `:clamp` (default) — value clamped to the nearest endpoint, a warning is
  emitted.
* `:line` — value linearly extrapolated using the slope of the nearest
  segment (equivalent to `Interpolations.LinearInterpolation(...; extrapolation_bc = Line())`).
  No warning. The same scheme is applied to the σ column (with `abs(...)` to keep
  uncertainties non-negative).

The file has either:

* **2 columns** `x, v` — the uncertainty on the value is taken as
  `err_pct * v` (default 10%); or
* **3 columns** `x, v, σ_v` — the third column carries the per-point 1-σ
  uncertainty and is interpolated alongside `v`.

Lines starting with `#` are treated as comments and ignored; leading non-numeric
rows are auto-detected and skipped. Use `skipstart = N` to force `N` physical
lines to be discarded from the top of the file before parsing.

Errors propagate by the standard division rule
`σ(y/v) = (y/v)·sqrt((σ_y/y)² + (σ_v/v)²)`, applied to both the peak yields and to
`tic`. Rows whose interpolated value is non-positive are left unchanged with a
warning.

(`normalize_flux` is a deprecated alias of this function.)
"""
function normalize_external(yc::YieldCurve, external_file::AbstractString;
                            err_pct::Real       = 0.10,
                            skipstart::Int      = 0,
                            extrapolate::Symbol = :clamp)
    extrapolate ∈ (:clamp, :line) ||
        error("normalize_external: extrapolate must be :clamp or :line (got :$extrapolate)")
    xf, ff, σf = _read_flux(external_file; skipstart = skipstart,
                            flux_err_pct = err_pct)
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
            @warn "normalize_external: x=$(yc.x[i]) outside reference range " *
                  "[$(xf[1]), $(xf[end])]; clamped to nearest" value = φ
        end
        if !(φ > 0)
            @warn "normalize_external: non-positive reference value at x=$(yc.x[i]); " *
                  "skipping division" value = φ
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
    md["normalize_external"]         = String(external_file)
    md["normalize_external_err_pct"] = Float64(err_pct)
    md["normalize_external_extrap"]  = String(extrapolate)
    return YieldCurve(copy(yc.x), yc.xlabel,
                      Y, Y_err, tic, tic_err,
                      copy(yc.found_mz), copy(yc.labels), copy(yc.windows),
                      copy(yc.files), md)
end

"""
    normalize_flux(yc, flux_file; flux_err_pct = 0.10, kw...)
Deprecated alias for [`normalize_external`](@ref) (the operation is general — the
external reference need not be a flux). `flux_err_pct` maps to `err_pct`.
"""
function normalize_flux(yc::YieldCurve, flux_file::AbstractString;
                        flux_err_pct::Real  = 0.10,
                        skipstart::Int      = 0,
                        extrapolate::Symbol = :clamp)
    Base.depwarn("`normalize_flux` is deprecated, use `normalize_external` " *
                 "(with `err_pct` instead of `flux_err_pct`).", :normalize_flux)
    return normalize_external(yc, flux_file; err_pct = flux_err_pct,
                              skipstart = skipstart, extrapolate = extrapolate)
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
