"""
Plotting module for MScontainer data type (MSscans and IonCurrent).
```julia-repl
julia> plot(scans[1])
julia> plot(chr)
julia> plot(scans[1], annot = :auto)                    # label the top peaks
julia> plot(scans[1], sequence = "PEPTIDE", ions = (:b, :y))
julia> plot(chr, annot = :auto)                         # label chromatographic peaks
```
"""
module plots

using Plots, RecipesBase   # used for plotting
using Printf               # runtime m/z formatting for annotation labels
using LinearAlgebra: diag, checksquare   # covariance-map marginal + square check
using Statistics: median                 # robust (median-based) colour limits


using MassJ:MScontainer
using MassJ:IonCurrent
using MassJ:maxic
using MassJ:MSscans
using MassJ:YieldCurve
using MassJ:stdev, sem
using MassJ:num2pnt
using MassJ:FragmentIon, fragment_ions
using MassJ:AbstractPeak, Peak, TargetPeak
using MassJ:ChromPeak, chrom_peaks
using MassJ:cut_autocorrelation

"""
    normalisation(ms::MassJ.MScontainer)
Normalization function for plotting mass spectra in relative intensity.
"""
function normalisation(ms::MScontainer)
    factor = 100. / ms.basePeakIntensity
    return ms.int .* factor
end

"""
    normalisation(cr::MassJ.IonCurrent)
Normalization function for plotting ion-current traces in relative intensity.
"""
function normalisation(cr::IonCurrent)
    factor = 100. / maxic(cr)
    return cr.ic .* factor
end

"""
    scaling(cr::MassJ.IonCurrent)
Abscissa values for display, returned as-is. Retention time is taken to be in
minutes — the convention the mzML loader and the v2 pipeline use — so it is not
rescaled here; drift time and compensation voltage are likewise unscaled.
"""
function scaling(cr::IonCurrent)
    return cr.x
end


# ----------------------------------------------------------------------------
# Annotation layer (round 4, phase 1)
# ----------------------------------------------------------------------------
#
# The analytical engines (`fragment_ions`, `chrom_peaks`, the peak descriptors,
# …) all return results keyed to an m/z (or a trace abscissa). The helpers below
# turn any such source into a uniform list of *candidate* labels, match each to
# the nearest data point, resolve overlaps, and hand the result to the recipes
# as plot annotations plus a (plotly) hover layer.
#
# Labels are drawn vertically (rotation = 90°). That reduces overlap to a 1-D
# minimum-gap problem on the abscissa, which a `RecipesBase` recipe *can* solve
# in data coordinates even though it cannot measure rendered text extent.
#
# A label is carried internally as a NamedTuple `(mz, text, charge, group, idx)`
# (the candidate produced by an adapter) and, once matched to a data point, as
# `(x, mz, y, text, charge, group)`. `text` may be empty — then the m/z itself is
# shown; `charge` is 0 when unknown; `group` drives `color_by`.

# Runtime printf of a single number with a user format string (default "%.2f").
_fmtnum(x::Real, fmt::AbstractString) = Printf.format(Printf.Format(fmt), x)

# Canonical fragment label without the charge suffix (charge is carried
# separately and rendered by `_disp` / `show_charge`).
function _frag_base(fi::FragmentIon)
    s = string(fi.series) * string(fi.position)
    fi.isomer == 1 && (s *= "′")
    fi.hshift != 0 && (s *= fi.hshift > 0 ? "+$(fi.hshift)" : string(fi.hshift))
    return s
end

# --- adapters: source -> Vector of candidate NamedTuples -------------------
# Each candidate is (mz, text, charge, group, idx). `idx > 0` means the source
# already points at a specific data index (the `:auto` case); `idx == 0` means
# "match to the nearest peak within tolerance".

_spectrum_candidates(src::Symbol, ms::MSscans) =
    (src === :auto || src === :peaks) ?
        [(mz = ms.mz[i], text = "", charge = 0, group = :peak, idx = i)
         for i in eachindex(ms.mz)] :
        error("unknown annotate symbol :$src (use :auto, a vector, or mz => text pairs)")

_spectrum_candidates(fi::FragmentIon, ms::MSscans) = _spectrum_candidates([fi], ms)

_spectrum_candidates(src::AbstractVector{<:FragmentIon}, ms::MSscans) =
    [(mz = fi.mz, text = _frag_base(fi), charge = fi.charge, group = fi.series, idx = 0)
     for fi in src]

_spectrum_candidates(p::AbstractPeak, ms::MSscans) = _spectrum_candidates([p], ms)

_spectrum_candidates(src::AbstractVector{<:AbstractPeak}, ms::MSscans) =
    [_peak_candidate(p) for p in src]

_peak_candidate(p::Peak) =
    (mz = 0.5 * (p.mz1 + p.mz2), text = p.label, charge = 0, group = :peak, idx = 0)
_peak_candidate(p::TargetPeak) =
    (mz = p.mz, text = p.label, charge = 0, group = :peak, idx = 0)

_spectrum_candidates(src::AbstractVector{<:Pair}, ms::MSscans) =
    [(mz = Float64(first(p)), text = string(last(p)), charge = 0, group = :peak, idx = 0)
     for p in src]

_spectrum_candidates(src::AbstractVector{<:Real}, ms::MSscans) =
    [(mz = Float64(m), text = "", charge = 0, group = :peak, idx = 0) for m in src]

# --- matching + de-clutter --------------------------------------------------

# Match every candidate to a data point, drop those absent from this spectrum,
# then greedily keep the most intense labels subject to a minimum m/z gap.
function _place_spectrum(ms::MSscans, cands, factor::Real;
                         tol::Real, ppm, declutter, nlabels::Integer)
    placed = NamedTuple[]
    for c in cands
        if c.idx > 0
            i = c.idx
        else
            i = num2pnt(ms.mz, c.mz)
            t = ppm === nothing ? Float64(tol) : c.mz * ppm * 1e-6
            abs(ms.mz[i] - c.mz) <= t || continue
        end
        push!(placed, (x = ms.mz[i], mz = ms.mz[i], y = ms.int[i] * factor,
                       text = c.text, charge = c.charge, group = c.group))
    end
    isempty(placed) && return placed
    sort!(placed, by = p -> p.y, rev = true)
    span = maximum(ms.mz) - minimum(ms.mz)
    return _declutter(placed, declutter, nlabels, span)
end

# Greedy 1-D suppression on the abscissa (`.x`). `declutter` is `:suppress`/`true`
# (default gap = span/40), `:none`/`false` (cap only), or an explicit gap.
function _declutter(placed, declutter, nlabels::Integer, span::Real)
    min_gap = (declutter === :none || declutter === false) ? 0.0 :
              declutter isa Real                            ? Float64(declutter) :
              span / 40
    out = empty(placed)
    for p in placed
        length(out) >= nlabels && break
        if min_gap <= 0 || all(o -> abs(o.x - p.x) >= min_gap, out)
            push!(out, p)
        end
    end
    return out
end

# Final display string for a placed spectrum label.
function _disp(p; show_mz::Bool, show_charge::Bool, fmt::AbstractString)
    s = p.text
    if show_charge && p.charge != 0 && abs(p.charge) != 1
        glyph = string(abs(p.charge), p.charge < 0 ? '-' : '+')
        s = isempty(s) ? glyph : s * " " * glyph
    end
    if show_mz || isempty(s)
        m = _fmtnum(p.mz, fmt)
        s = isempty(s) ? m : s * " " * m
    end
    return s
end

# Hover tooltip (plotly): always informative regardless of the on-plot text.
function _hover(p; fmt::AbstractString)
    h = "m/z " * _fmtnum(p.mz, fmt)
    p.charge != 0 && (h *= "  z=" * string(p.charge))
    isempty(p.text) || (h *= "  " * p.text)
    return h
end

# Per-label colours for `color_by` (:none -> all black; otherwise by group).
const _PALETTE = [:black, :red, :blue, :green, :purple, :orange, :brown, :magenta, :teal]
function _label_colors(placed, color_by::Symbol)
    color_by === :none && return fill(:black, length(placed))
    groups = unique(p.group for p in placed)
    cmap = Dict(g => _PALETTE[mod1(i, length(_PALETTE))] for (i, g) in enumerate(groups))
    return [cmap[p.group] for p in placed]
end

# Vertical text object placed just above a peak apex.
_textobj(s::AbstractString, c) =
    text(s; pointsize = 7, rotation = 90, halign = :left, valign = :bottom, color = c)


"""
    g(ms::MassJ.MSscans; method = :relative, band = :none, annot = nothing, …)

Plot a mass spectrum. Relative intensity by default; `method = :absolute` keeps
raw intensities. Set `band = :sem` / `:std` to draw a ±1-σ uncertainty ribbon on
an averaged spectrum (ignored when no replicate statistics are present).

Annotation (round 4): `annot` accepts `:auto` (the spectrum's most intense
peaks), a `Vector{FragmentIon}`, a `Vector{<:AbstractPeak}`, a bare m/z vector,
or `mz => "text"` pairs; pass `sequence = "PEPTIDE"` (with optional `ions` /
`charges`) as a shortcut for fragment-ion labels. Labels are drawn vertically
above the matched apex; overlaps are removed greedily (`declutter`). Keywords:
`nlabels` (max labels), `tol` / `ppm` (apex-match tolerance), `show_mz`,
`show_charge`, `annot_fmt` (m/z format), `declutter` (`:suppress` / `:none` / a
gap in m/z), `color_by` (`:none` / `:group`). `annot = nothing` (default)
reproduces the un-annotated plot exactly. Centroid spectra default to a
stem (`:sticks`) series. (`annot` rather than `annotate`, which Plots reserves
as an alias of its own `annotations` attribute.)
"""
@recipe function g(ms::MSscans; method = :relative, band = :none,
                   annot = nothing, sequence = nothing,
                   ions = (:b, :y), charges = 1:1,
                   nlabels = 10, tol = 0.5, ppm = nothing,
                   show_mz = false, show_charge = true, annot_fmt = "%.2f",
                   declutter = :suppress, color_by = :none)
    seriestype  --> (ms.spectrumType === :centroid ? :sticks : :path)
    seriescolor --> :red
    label       --> ""
    xlabel      --> "m/z"
    # base-peak height; fall back to the data max when the file/scan did not record
    # basePeakIntensity (e.g. mzML without the cvParam) so normalisation never divides by 0.
    base = ms.basePeakIntensity > 0 ? ms.basePeakIntensity :
           (isempty(ms.int) ? 1.0 : maximum(ms.int))
    if method == :relative
        factor = 100. / base
        ymax   = 100.0
        ylabel --> "Intensity (%)"
    else  # :absolute
        factor = 1.0
        ymax   = base
        ylabel --> "Intensity (a.u.)"
    end
    if (band === :std || band === :sem) && !isempty(ms.s)
        e = band === :sem ? sem(ms) : stdev(ms)
        ribbon    --> e .* factor
        fillalpha --> 0.15
    end

    # resolve the annotation source (sequence shortcut wins)
    src = sequence !== nothing ?
          fragment_ions(String(sequence); ions = ions, charges = charges) : annot
    if src !== nothing
        cands  = _spectrum_candidates(src, ms)
        placed = _place_spectrum(ms, cands, factor;
                                 tol = tol, ppm = ppm,
                                 declutter = declutter, nlabels = nlabels)
        if !isempty(placed)
            ylims --> (0.0, ymax * 1.3)
            cols   = _label_colors(placed, color_by)
            off    = 0.015 * ymax
            labels = [_textobj(_disp(p; show_mz = show_mz, show_charge = show_charge,
                                     fmt = annot_fmt), cols[i])
                      for (i, p) in enumerate(placed)]
            htxt   = [_hover(p; fmt = annot_fmt) for p in placed]
            @series begin
                seriestype         := :scatter
                primary            := false
                markersize         := 1
                markeralpha        := 0.0
                markerstrokewidth  := 0
                series_annotations := labels
                hover              := htxt
                [p.mz for p in placed], [p.y + off for p in placed]
            end
        end
    end

    ms.mz, ms.int .* factor
end


# --- IonCurrent annotation --------------------------------------------------

_trace_candidates(src::Symbol, cr::IonCurrent) =
    (src === :auto || src === :peaks) ? _trace_candidates(chrom_peaks(cr), cr) :
        error("unknown annotate symbol :$src (use :auto, a Vector{ChromPeak}, reals, or pairs)")

_trace_candidates(src::AbstractVector{<:ChromPeak}, ::IonCurrent) =
    [(x = p.apex, text = "") for p in src]
_trace_candidates(src::AbstractVector{<:Pair}, ::IonCurrent) =
    [(x = Float64(first(p)), text = string(last(p))) for p in src]
_trace_candidates(src::AbstractVector{<:Real}, ::IonCurrent) =
    [(x = Float64(v), text = "") for v in src]

function _place_trace(cr::IonCurrent, sx, y, cands;
                      declutter, nlabels::Integer, fmt::AbstractString)
    placed = NamedTuple[]
    for c in cands
        i  = num2pnt(cr.x, c.x)
        xx = sx[i]
        txt = isempty(c.text) ? _fmtnum(xx, fmt) : c.text
        push!(placed, (x = xx, y = y[i], text = txt))
    end
    isempty(placed) && return placed
    sort!(placed, by = p -> p.y, rev = true)
    span = isempty(sx) ? 0.0 : (maximum(sx) - minimum(sx))
    return _declutter(placed, declutter, nlabels, span)
end


"""
    h(cr::MassJ.IonCurrent; method = :relative, annot = nothing, …)

Plot an ion-current trace (chromatogram / mobilogram / ionogram). Relative
intensity by default; `method = :absolute` keeps raw intensities. Set
`annot = :auto` to label detected peaks ([`chrom_peaks`](@ref)), or pass a
`Vector{ChromPeak}`, a bare abscissa vector, or `x => "text"` pairs. Labels are
the apex abscissa unless given; same `nlabels` / `declutter` / `annot_fmt`
controls as the spectrum recipe. `annot = nothing` (default) is the plain
trace.
"""
@recipe function h(cr::IonCurrent; method = :relative,
                   annot = nothing, nlabels = 10,
                   annot_fmt = "%.2f", declutter = :suppress)
    seriestype  --> :path
    seriescolor --> :blue
    fillrange   --> 0
    fillalpha   --> 0.3
    label       --> ""
    xlabel      --> (cr.axis === :rt ? "time (mins)" :
                     cr.axis === :drift ? "drift time" : "compensation voltage")
    if method == :relative
        y = normalisation(cr)
        ymax = 100.0
        ylabel --> "Intensity (%)"
    elseif method == :absolute
        y = cr.ic
        ymax = isempty(cr.ic) ? 1.0 : maximum(cr.ic)
        ylabel --> "Intensity (a.u.)"
    end
    sx = scaling(cr)

    if annot !== nothing
        cands  = _trace_candidates(annot, cr)
        placed = _place_trace(cr, sx, y, cands;
                              declutter = declutter, nlabels = nlabels, fmt = annot_fmt)
        if !isempty(placed)
            ylims --> (0.0, ymax * 1.3)
            off    = 0.015 * ymax
            labels = [_textobj(p.text, :black) for p in placed]
            htxt   = [p.text for p in placed]
            @series begin
                seriestype         := :scatter
                primary            := false
                markersize         := 1
                markeralpha        := 0.0
                markerstrokewidth  := 0
                fillrange          := nothing
                series_annotations := labels
                hover              := htxt
                [p.x for p in placed], [p.y + off for p in placed]
            end
        end
    end

    sx, y
end


"""
    k(yc::MassJ.YieldCurve)
Plot a YieldCurve: one line per peak, x-axis from `yc.x`, legend from `yc.labels`,
x-label from `yc.xlabel`. When `yc.yields_err` carries finite uncertainties,
1-σ ribbons are drawn around each line with `fillalpha = 0.15`.
"""
@recipe function k(yc::YieldCurve)
    seriestype --> :path
    xlabel     --> yc.xlabel
    ylabel     --> "yield (a.u.)"
    label      --> reshape(yc.labels, 1, :)
    if any(isfinite, yc.yields_err)
        ribbon    --> replace(yc.yields_err, NaN => 0.0)
        fillalpha --> 0.15
    end
    yc.x, yc.yields
end


# --- peak-window overlay ----------------------------------------------------

# Flat [lo1, hi1, lo2, hi2, …] window bounds for a :vspan series.
_peak_bounds(p::Peak) = [p.mz1, p.mz2]
function _peak_bounds(p::TargetPeak)
    b = Float64[]
    for m in p.mzs
        push!(b, m - p.tol, m + p.tol)
    end
    return b
end

"""
    pk(p::MassJ.AbstractPeak)

Overlay a peak's integration window(s) on an existing spectrum plot as shaded
vertical spans, so you can see where a [`Peak`](@ref) / [`TargetPeak`](@ref)
integrates. For a multi-target `TargetPeak` (isotope cluster) every sub-window is
drawn. Use `plot!` after plotting a spectrum:

```julia
plot(avg)
plot!(TargetPeak("Nd(NO3)4", "Precursor"; charge = -1, tol = 0.2))
```

The span colour/transparency/label follow the usual Plots keywords
(`seriescolor`, `seriesalpha`, `label`).
"""
@recipe function pk(p::AbstractPeak)
    seriestype  --> :vspan
    seriesalpha --> 0.2
    seriescolor --> :orange
    label       --> p.label
    _peak_bounds(p)
end

# Palette for overlaying several peak windows (avoids red — the spectrum colour).
const _SPAN_PALETTE = [:orange, :green, :blue, :purple, :brown, :magenta, :teal, :olive, :cyan]

"""
    pkv(ps::AbstractVector{<:MassJ.AbstractPeak})

Overlay a whole list of peaks at once, each window shaded with a colour cycled
from a palette and labelled by `peak.label`:

```julia
plot(avg)
plot!(peaks)          # all TargetPeak / Peak windows in one call
```

Equivalent to `plot!`-ing each peak in turn; see `pk` for a single peak.
"""
@recipe function pkv(ps::AbstractVector{<:AbstractPeak})
    seriesalpha --> 0.2
    for (i, p) in enumerate(ps)
        @series begin
            seriestype  := :vspan
            seriescolor --> _SPAN_PALETTE[mod1(i, length(_SPAN_PALETTE))]
            label       := p.label
            _peak_bounds(p)
        end
    end
end


export association_map, association_map_marginal, significance_map

# Colour scheme + limits for an association map. Colour: a diverging scale (RdBu) for a
# signed map (covariance, correlation), a sequential scale (viridis) from zero for a
# non-negative one (conditional mutual information). Limits: unless `clims` is given
# explicitly, the amplitude is *robust* — `climscale` times the larger of the median
# positive / median negative magnitude of the off-centre values — so that a handful of
# strong pairs saturate (staying fully visible) instead of stretching the range and
# washing everything else out (the failure mode of a max-based scale on a sparse map). For
# a signed map the limits are kept **symmetric** `(-amp, amp)` so that zero maps to the
# neutral (white) colour. Larger `climscale` widens the range (less saturation, more
# dynamic range); smaller narrows it (more sensitive to weak associations). Non-finite
# (masked / diagonal-cut) cells are ignored throughout.
function _map_colorscale(M; color = nothing, clims = nothing, climscale::Real = 6.0)
    vals = filter(isfinite, vec(M))
    signed = !isempty(vals) && minimum(vals) < 0
    color === nothing && (color = signed ? :RdBu : :viridis)
    if clims === nothing
        if isempty(vals)
            clims = signed ? (-1.0, 1.0) : (0.0, 1.0)
        elseif signed
            pos = filter(>(0), vals); neg = filter(<(0), vals)
            mp  = isempty(pos) ? 0.0 : median(pos)
            mn  = isempty(neg) ? 0.0 : abs(median(neg))
            amp = climscale * max(mp, mn)
            amp > 0 || (amp = maximum(abs, vals); amp = amp > 0 ? amp : 1.0)   # degenerate
            clims = (-amp, amp)                             # symmetric → 0 is white
        else
            pos = filter(>(0), vals)
            hi = isempty(pos) ? maximum(vals) : climscale * median(pos)
            clims = (0.0, hi > 0 ? float(hi) : 1.0)
        end
    end
    return color, clims
end

# Copy `A` to a Float matrix, blank the autocorrelation diagonal band (to `NaN`) when
# `cut`, and blank cells where `mask` is false (to `NaN`) — the shared preparation for
# the map plotters. `NaN` cells render transparent and are ignored by `_map_colorscale`.
function _prepare_map(A, n; cut, halfwidth, mask)
    mask === nothing || Base.size(mask) == (n, n) ||
        error("association map: mask must be $(n)×$(n) to match the map.")
    M = cut ? cut_autocorrelation(float.(A); halfwidth = halfwidth, value = NaN) : float.(A)
    if mask !== nothing
        @inbounds for i in eachindex(M, mask)
            mask[i] || (M[i] = NaN)
        end
    end
    return M
end

"""
    association_map(A, mz = 1:size(A,1); cut = true, halfwidth = 1, mask = nothing,
                    climscale = 6, color = nothing, clims = nothing, kwargs...)

Heatmap of a square association map `A` — a covariance (`cov`),
correlation ([`MassJ.pcor`](@ref)), or conditional-mutual-information
([`MassJ.cmi`](@ref)) matrix — with its columns indexed by `mz`. A signed map gets a
diverging colour scale (positive = same precursor, negative = different); a non-negative
map (CMI) gets a sequential scale from zero. `cut` suppresses the autocorrelation diagonal
band (half-width `halfwidth`) so cross-peaks are visible.

By default the colour limits are **robust and symmetric**: the amplitude is `climscale`
times the larger of the median positive / median negative association, applied as
`(-amp, amp)` so zero maps to the neutral (white) colour. A few strong cross-peaks then
saturate (staying visible) instead of stretching the scale and washing the map out. Raise
`climscale` for more dynamic range, lower it to bring up weak associations; or pass
`clims = (lo, hi)` to set the range by hand (which overrides `climscale`). Returns a `Plots.Plot`; overlay picked pairs with
`scatter!` on the `MapFeature` `mz1`/`mz2`. For the framed figure with marginal spectra see
[`association_map_marginal`](@ref).

Pass `mask` — a `Bool` matrix the same size as `A`, e.g. the FDR gate from
[`MassJ.fdr_adjust`](@ref) or a [`MassJ.correspondence`](@ref) result's `significant` field — to keep
only the edges where `mask` is true; the rest render transparent. This turns a raw score
map into a *significance-gated* map. To plot the p-values themselves see
[`significance_map`](@ref).

```julia
using MassJ.plots
association_map(MassJ.cov(X), mz)                      # or a pcor / cmi map
jk = jackknife(pcor, am.matrix, am.tic)
association_map(jk.estimate, am.mz; mask = fdr_adjust(jk.pvalue).significant)   # gated map
```
"""
function association_map(A::AbstractMatrix{<:Real},
                         mz::AbstractVector{<:Real} = 1:size(A, 1);
                         cut::Bool = true, halfwidth::Integer = 1,
                         mask::Union{Nothing,AbstractMatrix{Bool}} = nothing,
                         climscale::Real = 6, color = nothing, clims = nothing, kwargs...)
    n = checksquare(A)
    length(mz) == n || error("association_map: length(mz) must equal size(A,1) = $n.")
    M = _prepare_map(A, n; cut = cut, halfwidth = halfwidth, mask = mask)
    c, cl = _map_colorscale(M; color = color, clims = clims, climscale = climscale)
    return heatmap(mz, mz, M; color = c, clims = cl, xguide = "m/z", yguide = "m/z",
                   colorbar_title = "association", kwargs...)
end


"""
    association_map_marginal(A, mz, intensity = sqrt.(max.(diag(A), 0)); cut = true,
                             halfwidth = 1, mask = nothing, climscale = 6,
                             color = nothing, clims = nothing, size = (720, 720), kwargs...)

The classic covariance-mapping figure: the square association map `A` as a central
heatmap (coloured as in [`association_map`](@ref)), with the 1-D marginal spectrum
`intensity` (length `size(A,1)`) drawn *above* it (sharing the m/z x-axis) and to its
*left* (sharing the m/z y-axis), so a cross-peak lines up with the two ions that
produced it. `intensity` defaults to the per-bin standard deviation (√ of `A`'s
diagonal), which is meaningful for a *covariance* map; for a correlation or CMI map
pass the averaged or summed series spectrum instead. `cut` suppresses the
autocorrelation band (half-width `halfwidth`); `mask` gates the heatmap to significant
edges, and `climscale`/`clims` set the (robust, median-based by default) colour range —
all exactly as in [`association_map`](@ref). Returns a `Plots.Plot` and works with
any backend (`gr()` for print, `plotly()` for interactive).

```julia
using MassJ.plots
am = MassJ.abundance_matrix(scans; binsize = 1.0)
A  = MassJ.cov(am.matrix ./ max.(am.tic, 1))   # TIC-normalised
association_map_marginal(A, am.mz, vec(sum(am.matrix, dims = 1)))
```
"""
function association_map_marginal(A::AbstractMatrix{<:Real}, mz::AbstractVector{<:Real},
                                  intensity::AbstractVector{<:Real} = sqrt.(max.(diag(A), 0.0));
                                  cut::Bool = true, halfwidth::Integer = 1,
                                  mask::Union{Nothing,AbstractMatrix{Bool}} = nothing,
                                  climscale::Real = 6, color = nothing, clims = nothing,
                                  size = (720, 720), kwargs...)
    n = checksquare(A)
    (length(mz) == n && length(intensity) == n) ||
        error("association_map_marginal: mz and intensity must have length size(A,1) = $n.")
    M = _prepare_map(A, n; cut = cut, halfwidth = halfwidth, mask = mask)
    c, cl = _map_colorscale(M; color = color, clims = clims, climscale = climscale)
    xl = (float(first(mz)), float(last(mz)))

    top = plot(mz, intensity; xlims = xl, legend = false, grid = false,
               framestyle = :box, ticks = false, linecolor = :steelblue, fillrange = 0,
               fillalpha = 0.15, fillcolor = :steelblue)
    left = plot(intensity, mz; ylims = xl, xflip = true, legend = false, grid = false,
                framestyle = :box, ticks = false, linecolor = :steelblue, fillrange = 0,
                fillalpha = 0.15, fillcolor = :steelblue)
    hm = heatmap(mz, mz, M; color = c, clims = cl, xlims = xl, ylims = xl,
                 xguide = "m/z", yguide = "m/z", colorbar_title = "association")

    lay = @layout [ _         t{0.80w,0.20h}
                    s{0.20w}  h ]
    return plot(top, left, hm; layout = lay, size = size, legend = false, kwargs...)
end


"""
    significance_map(P, mz = 1:size(P,1); q = nothing, transform = :neglog10,
                     cut = true, halfwidth = 1, color = nothing, clims = nothing, kwargs...)

Heatmap of a *p-value* map `P` — the `pvalue` field returned by [`MassJ.jackknife`](@ref) or
[`MassJ.permutation`](@ref) — transformed so that strong evidence reads bright. With
`transform = :neglog10` (default) each cell is shown as ``-\\log_{10} p`` on a sequential
scale (large = significant); `transform = :none` plots `P` itself. When `q` is given,
cells with ``p \\ge q`` are blanked (rendered transparent), leaving only the edges that
survive an [`MassJ.fdr_adjust`](@ref)-style threshold. `cut` suppresses the (untested)
autocorrelation diagonal. Columns are indexed by `mz`. Returns a `Plots.Plot`.

To instead colour a gated map by the *score* (rather than the p-value), pass the gate as
`mask` to [`association_map`](@ref).

```julia
using MassJ.plots
jk = jackknife(pcor, am.matrix, am.tic)
significance_map(jk.pvalue, am.mz; q = 0.05)
```
"""
function significance_map(P::AbstractMatrix{<:Real},
                          mz::AbstractVector{<:Real} = 1:size(P, 1);
                          q::Union{Nothing,Real} = nothing, transform::Symbol = :neglog10,
                          cut::Bool = true, halfwidth::Integer = 1,
                          color = nothing, clims = nothing, kwargs...)
    n = checksquare(P)
    length(mz) == n || error("significance_map: length(mz) must equal size(P,1) = $n.")
    all(p -> isnan(p) || (0 <= p <= 1), P) ||
        error("significance_map: P must contain p-values in [0, 1].")
    V = transform === :neglog10 ? -log10.(clamp.(float.(P), 1e-300, 1.0)) :
        transform === :none     ? float.(P) :
        error("significance_map: transform must be :neglog10 or :none.")
    if q !== nothing
        @inbounds for i in eachindex(V, P)
            (isnan(P[i]) || P[i] >= q) && (V[i] = NaN)
        end
    end
    Vc = cut ? cut_autocorrelation(V; halfwidth = halfwidth, value = NaN) : V
    c, cl = _map_colorscale(Vc; color = color, clims = clims)
    ctitle = transform === :neglog10 ? "-log10 p" : "p-value"
    return heatmap(mz, mz, Vc; color = c, clims = cl, xguide = "m/z", yguide = "m/z",
                   colorbar_title = ctitle, kwargs...)
end


end # submodule
