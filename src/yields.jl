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
       read_peaklist, drop_yields, keep_yields, sum_yields, group_yields, with_yields,
       drop_peaks,
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


# --- Peak selection, summing, grouping, and derived columns -----------------
#
# The yields matrix is plain data: the intended workflow for arbitrary arithmetic
# is "select the columns, compute in Julia, put the result back" (`with_yields`,
# with Measurements.jl propagating the errors). The functions here cover the
# bookkeeping that is easy to get wrong when done by hand: selecting columns by
# label pattern, and merging co-named peaks (an ion series, a neutral loss) with a
# consistent error rule. All return a new, immutable `YieldCurve`; `tic`/`tic_err`
# are never touched (they keep totalling the *original* peak set, which is what
# `normalize_tic` should divide by); every operation is logged in
# `metadata["yield_ops"]`.

# A *peak selector* names columns by label: a single label, a collection of labels,
# a `Regex` (labels that match), or a predicate `label -> Bool`. `_peak_indices`
# resolves it to column indices. With `strict = true` a label absent from
# `yc.labels` is an error (a silent miss would corrupt a sum); `drop_yields` is
# lenient, as it always was.
_peak_indices(yc::YieldCurve, sel::AbstractString; strict::Bool = true) =
    _peak_indices(yc, (sel,); strict = strict)
_peak_indices(yc::YieldCurve, sel::Regex; strict::Bool = true) =
    findall(l -> occursin(sel, l), yc.labels)
_peak_indices(yc::YieldCurve, sel::Function; strict::Bool = true) =
    findall(l -> sel(l)::Bool, yc.labels)
function _peak_indices(yc::YieldCurve, sel; strict::Bool = true)
    idx = Int[]
    for s in sel
        s isa AbstractString ||
            error("peak selector: expected a label, a collection of labels, a Regex, or a " *
                  "predicate (got an element of type $(typeof(s))).")
        j = findfirst(==(String(s)), yc.labels)
        if j === nothing
            strict && error("peak selector: no peak labelled \"$s\" " *
                            "(labels: $(join(yc.labels, ", "))).")
        elseif !(j in idx)
            push!(idx, j)
        end
    end
    return idx
end

# Append one record to the `"yield_ops"` provenance list of a metadata dict (the
# list is copied first so the source curve's own list is never mutated).
function _note_peak_op!(md::Dict{String,Any}, op::AbstractString, label, sources)
    ops = copy(get(md, "yield_ops", Any[]))
    push!(ops, (op = String(op),
                label = label === nothing ? nothing : String(label),
                sources = sources === nothing ? nothing : String.(collect(sources))))
    md["yield_ops"] = ops
    return md
end

# Build a new YieldCurve keeping the columns whose indices are in `keep` (in that
# order) — the column analogue of `_take_rows`. `x`, `tic`, `tic_err`, `files` are
# untouched.
function _take_cols(yc::YieldCurve, keep::AbstractVector{<:Integer},
                    op::AbstractString, sources)
    md = _note_peak_op!(copy(yc.metadata), op, nothing, sources)
    return YieldCurve(copy(yc.x), yc.xlabel,
                      yc.yields[:, keep], yc.yields_err[:, keep],
                      copy(yc.tic), copy(yc.tic_err),
                      yc.found_mz[:, keep], yc.labels[keep], yc.windows[keep],
                      copy(yc.files), md)
end

"""
    keep_yields(yc::YieldCurve, sel) -> YieldCurve

Return a new [`YieldCurve`](@ref) holding only the peaks selected by `sel`, in the
order they are selected — the inverse of [`drop_yields`](@ref). A **peak selector**
`sel` is any of:

* a label, `"y5"`, or a collection of labels, `["a3", "b3"]` (that order is kept);
* a `Regex` matched against the labels, `r"^y\\d"` (all y-ions), `r"-H2O\$"`;
* a predicate on the label, `l -> startswith(l, "b")`.

A label absent from `yc.labels` is an error (unlike `drop_yields`, which ignores
it). `tic`/`tic_err` are left unchanged, so they still total the original peak
set. The same selectors are accepted by [`drop_yields`](@ref),
[`sum_yields`](@ref), and [`group_yields`](@ref); `yc["y5"]` returns one peak's
columns directly.

```julia
plot(keep_yields(yc, r"^y\\d"))                 # the y-ion series only
keep_yields(yc, ["precursor", "b3"])            # two peaks, in that order
```
"""
function keep_yields(yc::YieldCurve, sel)
    keep = _peak_indices(yc, sel; strict = true)
    return _take_cols(yc, keep, "keep", yc.labels[keep])
end

"""
    drop_yields(yc::YieldCurve, sel) -> YieldCurve

Return a new [`YieldCurve`](@ref) with the peaks selected by `sel` removed — a
label, a collection of labels, a `Regex`, or a predicate on the label (see
[`keep_yields`](@ref) for the selector forms). Labels not present in `yc.labels`
are silently ignored.

The `tic` field is left unchanged so it still reflects the totals over the
*original* peak set — useful when you want to keep the same TIC reference
after dropping a dominant peak (e.g. for plotting fragment yields without
the precursor swamping the axes).

```julia
plot(drop_yields(yc, "precursor"))
plot(drop_yields(yc, ["precursor", "solvent"]))
plot(drop_yields(yc, r"-H2O\$"))                # every water-loss channel
```
"""
function drop_yields(yc::YieldCurve, sel)
    drop = _peak_indices(yc, sel; strict = false)
    keep = setdiff(1:length(yc.labels), drop)
    return _take_cols(yc, keep, "drop", yc.labels[drop])
end

"""
    drop_peaks(yc::YieldCurve, sel)
Deprecated alias for [`drop_yields`](@ref). The column operations on a
[`YieldCurve`](@ref) are named `*_yields` (matching [`combine_yields`](@ref),
[`scale_yields`](@ref), [`trim_yields`](@ref)), leaving `*_peaks` to the functions
that build or detect [`Peak`](@ref)s — `fragment_peaks`, [`read_peaklist`](@ref),
[`chrom_peaks`](@ref).
"""
function drop_peaks(yc::YieldCurve, sel)
    Base.depwarn("`drop_peaks` is deprecated, use `drop_yields`.", :drop_peaks)
    return drop_yields(yc, sel)
end

"""
    yc[label] -> (y, err)

Index a [`YieldCurve`](@ref) by peak label: the yield column and its 1-σ column
as a named tuple, e.g. `y, err = yc["y5"]`. Use [`keep_yields`](@ref) to select
several peaks as a curve.
"""
function Base.getindex(yc::YieldCurve, label::AbstractString)
    j = findfirst(==(String(label)), yc.labels)
    j === nothing &&
        error("YieldCurve: no peak labelled \"$label\" (labels: $(join(yc.labels, ", "))).")
    return (y = yc.yields[:, j], err = yc.yields_err[:, j])
end

# Quadrature combination of 1-σ columns: NaN entries are skipped, and the result
# is NaN only when no member carries an error — the rule `tic_err` follows.
function _quadrature_err(E::AbstractMatrix, cols)
    n = size(E, 1)
    out = fill(NaN, n)
    @inbounds for i in 1:n
        acc = 0.0; any_err = false
        for j in cols
            e = E[i, j]
            if isfinite(e); acc += e * e; any_err = true; end
        end
        any_err && (out[i] = sqrt(acc))
    end
    return out
end

# Merge groups of columns into one column each. `groups` is a list of
# `name => column indices`. Every group column takes the position of its first
# member and, unless `keep_sources`, replaces its members; columns outside every
# group are untouched. Values add exactly (a yield is a mean of per-scan areas);
# errors add in quadrature — an independence assumption across the merged peaks.
# `found_mz` of a group column is NaN and its window is the envelope of the
# members'.
function _merge_cols(yc::YieldCurve, groups::AbstractVector, op::AbstractString;
                     keep_sources::Bool = false)
    n, m  = size(yc.yields)
    owner = zeros(Int, m)                       # column -> group index (0 = none)
    for (g, (name, cols)) in enumerate(groups)
        isempty(cols) && error("$(op)_yields: group \"$name\" selects no peak.")
        for j in cols
            owner[j] == 0 ||
                error("$(op)_yields: peak \"$(yc.labels[j])\" belongs to two groups " *
                      "(\"$(first(groups[owner[j]]))\" and \"$name\").")
            owner[j] = g
        end
    end
    Y = Vector{Float64}[]; E = Vector{Float64}[]; F = Vector{Float64}[]
    labels = String[]; windows = Tuple{Float64,Float64}[]
    emitted = falses(length(groups))
    for j in 1:m
        g = owner[j]
        if g != 0 && !emitted[g]
            name, cols = groups[g]
            push!(Y, vec(sum(view(yc.yields, :, cols), dims = 2)))
            push!(E, _quadrature_err(yc.yields_err, cols))
            push!(F, fill(NaN, n))
            push!(labels, String(name))
            push!(windows, (minimum(first, yc.windows[cols]), maximum(last, yc.windows[cols])))
            emitted[g] = true
        end
        if g == 0 || keep_sources
            push!(Y, yc.yields[:, j]); push!(E, yc.yields_err[:, j]); push!(F, yc.found_mz[:, j])
            push!(labels, yc.labels[j]); push!(windows, yc.windows[j])
        end
    end
    if !allunique(labels)
        dup = unique(Base.filter(l -> count(==(l), labels) > 1, labels))
        error("$(op)_yields: the result would carry duplicate labels $dup; choose another " *
              "group name / `label`, or drop the clashing peak first.")
    end
    md = copy(yc.metadata)
    for (name, cols) in groups
        _note_peak_op!(md, op, name, yc.labels[cols])
    end
    empty = Matrix{Float64}(undef, n, 0)
    return YieldCurve(copy(yc.x), yc.xlabel,
                      reduce(hcat, Y; init = empty), reduce(hcat, E; init = empty),
                      copy(yc.tic), copy(yc.tic_err),
                      reduce(hcat, F; init = empty), labels, windows, copy(yc.files), md)
end

"""
    sum_yields(yc::YieldCurve[, sel]; exclude = nothing, label = "sum", keep_sources = false)

Add peak columns of a [`YieldCurve`](@ref) into one column: every column when
`sel` is omitted, the peaks selected by `sel` (a label, a collection of labels, a
`Regex`, or a predicate — see [`keep_yields`](@ref)), or every column *except* those
selected by `exclude`. The summed column, named `label`, takes the position of its
first member and replaces the members unless `keep_sources = true`; the other
columns are untouched.

Values add exactly (a yield is a mean of per-scan areas, and the mean is linear).
The 1-σ column is the **quadrature** sum of the members' errors — the same rule
`tic_err` follows — which assumes the merged peaks fluctuate independently across
scans. For peaks that co-vary (isotopologues, fragments of one precursor) this is
only approximate; the rigorous route is a cluster [`TargetPeak`](@ref), which sums
the windows *inside each scan* before the error is taken. `tic`/`tic_err` are
left unchanged (they still total the original peak set), `found_mz` of the new
column is `NaN`, and its `windows` entry is the envelope of the members'. The
operation is recorded in `metadata["yield_ops"]`.

```julia
sum_yields(yc)                                          # everything → "sum"
sum_yields(yc, r"^a\\d"; label = "a-ions")              # the a-ion series
sum_yields(yc; exclude = "precursor", label = "fragments")
```
See [`group_yields`](@ref) for several sums at once and [`with_yields`](@ref) for
any other combination.
"""
function sum_yields(yc::YieldCurve, sel = nothing; exclude = nothing,
                   label::AbstractString = "sum", keep_sources::Bool = false)
    (sel !== nothing && exclude !== nothing) &&
        error("sum_yields: pass either a selector or `exclude`, not both.")
    m = length(yc.labels)
    cols = sel !== nothing     ? _peak_indices(yc, sel; strict = true) :
           exclude !== nothing ? setdiff(1:m, _peak_indices(yc, exclude; strict = true)) :
                                 collect(1:m)
    isempty(cols) && error("sum_yields: the selection is empty.")
    return _merge_cols(yc, [String(label) => cols], "sum"; keep_sources = keep_sources)
end

# Group name for one label under a function / Regex key (`nothing` = not grouped).
_group_name(f::Function, l::AbstractString) = (r = f(l); r === nothing ? nothing : string(r))
function _group_name(re::Regex, l::AbstractString)
    m = match(re, l)
    m === nothing && return nothing
    return isempty(m.captures) ? String(m.match) : String(something(m.captures[1], m.match))
end

# Function / Regex keys: one pass over the labels, groups in first-appearance order.
function _resolve_groups(yc::YieldCurve, key::Union{Function,Regex})
    names = String[]; cols = Vector{Int}[]
    for (j, l) in enumerate(yc.labels)
        name = _group_name(key, l)
        name === nothing && continue
        g = findfirst(==(name), names)
        g === nothing ? (push!(names, name); push!(cols, [j])) : push!(cols[g], j)
    end
    return [names[g] => cols[g] for g in eachindex(names)]
end

# Explicit `name => selector` pairs.
function _resolve_groups(yc::YieldCurve, pairs)
    out = Pair{String,Vector{Int}}[]
    for p in pairs
        p isa Pair ||
            error("group_yields: `key` must be a function, a Regex, or `name => selector` pairs.")
        push!(out, String(first(p)) => _peak_indices(yc, last(p); strict = true))
    end
    return out
end

"""
    group_yields(yc::YieldCurve, key; keep_sources = false)

Merge the peaks of a [`YieldCurve`](@ref) into groups named by a criterion on the
labels, producing one summed column per group (each computed as in
[`sum_yields`](@ref)). `key` is one of:

* a function `label -> name`, returning the group name (a `String` or `Symbol`),
  or `nothing` to leave that peak alone;
* a `Regex`: the first capture group names the group (the whole match when the
  pattern has no capture group); non-matching peaks are left alone;
* a vector of `name => selector` pairs (selectors as in [`keep_yields`](@ref)) for
  hand-named groups. Every group must select at least one peak, and a peak may
  belong to one group only.

Groups appear in order of first appearance (function / `Regex`) or as given
(pairs); each group column takes the position of its first member, and the members
are replaced unless `keep_sources = true`. Everything said in [`sum_yields`](@ref)
about errors, `tic`, `found_mz`, `windows`, and provenance applies per group.

With the labels produced by `fragment_peaks` (`"b3"`, `"y5(2+)"`, …) and
hand-labelled neutral losses such as `"b3-H2O"`:

```julia
group_yields(yc, r"^([abxy])\\d")                  # one column per ion series
group_yields(yc, r"-(\\w+)\$")                      # one column per neutral loss (H2O, NH3, …)
group_yields(yc, r"\\((\\d)\\+\\)")                   # by charge state
group_yields(yc, l -> occursin('(', l) ? "multiply charged" : nothing)
group_yields(yc, ["N-term" => r"^[abc]\\d", "C-term" => r"^[xyz]\\d"])
```
"""
function group_yields(yc::YieldCurve, key; keep_sources::Bool = false)
    groups = _resolve_groups(yc, key)
    isempty(groups) && error("group_yields: the key matched no peak.")
    return _merge_cols(yc, groups, "group"; keep_sources = keep_sources)
end

"""
    with_yields(yc::YieldCurve, labels, Y[, Yerr]; replace = false, sources = nothing,
               windows = nothing, found_mz = nothing) -> YieldCurve

Put columns computed *outside* the yields machinery back into a
[`YieldCurve`](@ref) — the way back after "take the matrix, compute in Julia".
`labels` is a label or a vector of labels and `Y` the matching yield vector or
`n × k` matrix (one column per label, rows aligned with `yc.x`); `Yerr` carries the
1-σ values in the same shape (`NaN` when omitted). The new columns are appended;
with `replace = true` a label that already exists is overwritten in place instead
(keeping its position, `windows`, and `found_mz`), which is how a column is
transformed in situ. Without `replace`, an existing label is an error.

`tic`/`tic_err` are left unchanged (they total the original peak set). The
operation is recorded in `metadata["yield_ops"]`; pass `sources` (a vector of
label vectors, one per new column — or a flat label vector for a single column) to
record what each column was derived from. `windows` (a vector of `(lo, hi)`
tuples) and `found_mz` (an `n × k` matrix) may be supplied for appended columns;
they default to `(NaN, NaN)` and `NaN`.

With `Measurements` loaded, `Y` may be a vector or matrix of `Measurement`s — the
value and uncertainty are split automatically, so any expression built on
[`measurements`](@ref)`(yc)` round-trips with its propagated 1-σ:

```julia
using Measurements
M = measurements(yc)                         # yield ± err, columns follow yc.labels
a, b = M[:, 1], M[:, 2]
yc2 = with_yields(yc, ["a − b", "a fraction"], hcat(a .- b, a ./ (a .+ b)))
yc3 = with_yields(yc, "y5", sqrt.(M[:, 3]); replace = true)   # in place
```
Without `Measurements`:
```julia
d = yc.yields[:, 1] .- yc.yields[:, 2]
σ = sqrt.(yc.yields_err[:, 1].^2 .+ yc.yields_err[:, 2].^2)
with_yields(yc, "a − b", d, σ; sources = ["a1", "b1"])
```
"""
function with_yields(yc::YieldCurve, labels, Y::AbstractVecOrMat{<:Real},
                    Yerr::Union{AbstractVecOrMat{<:Real},Nothing} = nothing;
                    replace::Bool = false, sources = nothing,
                    windows = nothing, found_mz = nothing)
    n    = length(yc.x)
    labs = labels isa AbstractString ? [String(labels)] : String.(collect(labels))
    k    = length(labs)
    Ym   = Y isa AbstractVector ? reshape(collect(Float64, Y), :, 1) : Matrix{Float64}(Y)
    Em   = Yerr === nothing         ? fill(NaN, size(Ym)) :
           Yerr isa AbstractVector  ? reshape(collect(Float64, Yerr), :, 1) :
                                      Matrix{Float64}(Yerr)
    size(Ym, 1) == n || error("with_yields: Y has $(size(Ym, 1)) rows, expected $n (one per x).")
    size(Ym, 2) == k || error("with_yields: Y has $(size(Ym, 2)) column(s) for $k label(s).")
    size(Em) == size(Ym) || error("with_yields: Yerr must have the same shape as Y.")
    allunique(labs) || error("with_yields: labels must be unique.")
    src = sources === nothing ? fill(nothing, k) :
          (k == 1 && all(s -> s isa AbstractString, sources)) ? Any[sources] :
          length(sources) == k ? collect(sources) :
          error("with_yields: `sources` must have one entry (a vector of labels) per new column.")
    Fm  = found_mz === nothing ? fill(NaN, n, k) :
          (size(found_mz) == (n, k) ? Matrix{Float64}(found_mz) :
           error("with_yields: `found_mz` must be an $n × $k matrix."))
    Wv  = windows === nothing ? fill((NaN, NaN), k) :
          (length(windows) == k ? [(Float64(w[1]), Float64(w[2])) for w in windows] :
           error("with_yields: `windows` must hold one (lo, hi) tuple per new column."))

    Yout = copy(yc.yields); Eout = copy(yc.yields_err); Fout = copy(yc.found_mz)
    lout = copy(yc.labels); wout = copy(yc.windows)
    md   = copy(yc.metadata)
    new  = Int[]
    for (c, l) in enumerate(labs)
        j = findfirst(==(l), yc.labels)
        if j !== nothing && replace
            Yout[:, j] = Ym[:, c]; Eout[:, j] = Em[:, c]
            _note_peak_op!(md, "replace", l, src[c])
        elseif j !== nothing
            error("with_yields: label \"$l\" already exists; pass `replace = true` to overwrite " *
                  "it in place, or drop it first.")
        else
            push!(new, c)
            _note_peak_op!(md, "with", l, src[c])
        end
    end
    if !isempty(new)
        Yout = hcat(Yout, Ym[:, new]); Eout = hcat(Eout, Em[:, new]); Fout = hcat(Fout, Fm[:, new])
        append!(lout, labs[new]); append!(wout, Wv[new])
    end
    return YieldCurve(copy(yc.x), yc.xlabel, Yout, Eout, copy(yc.tic), copy(yc.tic_err),
                      Fout, lout, wout, copy(yc.files), md)
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
# (`trim_yields` and `restrict_x` are the row analogue of `drop_yields`.)

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
