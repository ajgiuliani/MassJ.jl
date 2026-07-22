# Plotting
Plotting facilities are available as a submodule to the `MassJ` package. The [`MassJ.plots`](@ref) module relies on the [RecipesBase package](https://github.com/JuliaPlots/RecipesBase.jl), which allows writing recipes to plot users' data types. Recipes are provided for `MSscans`, `MSscans`, `IonCurrent`, and `YieldCurve`:

```julia
plot(scans[1], method = :relative)
plot(yc)                                   # YieldCurve: one line per peak
```

For mass spectra and chromatograms, plotting is made in relative intensities by default; this can be changed by setting `method = :absolute`.

For an averaged mass spectrum, a ±1-σ uncertainty ribbon can be drawn by passing `band = :sem` (standard error of the mean) or `band = :std` (sample standard deviation). The ribbon is scaled consistently in both `:relative` and `:absolute` modes and is silently ignored on a single scan that carries no replicate statistics. See [Uncertainties and errors](@ref) for the accessors behind it.

```julia
plot(avg; band = :sem)                     # ±1-σ (standard error of the mean)
plot(avg; band = :std, method = :absolute) # ±1-σ (sample standard deviation)
```

## Peak annotation

Spectra and ion-current traces can be annotated directly from the analytical
engines. Pass `annot` to label peaks; the labels are drawn vertically just above
the matched apex, and overlapping labels are removed greedily so the figure stays
readable. `annot = nothing` (the default) reproduces the plain plot.

```julia
plot(scan, annot = :auto)                        # the most intense peaks (m/z labels)
plot(scan, annot = [195.09 => "M+H", 217.07 => "M+Na"])  # explicit mz => text
plot(scan, annot = peaks)                         # a Vector of Peak / TargetPeak
plot(chrom, annot = :auto)                        # chromatographic peaks (chrom_peaks)
```

`annot` accepts `:auto` (top peaks of the spectrum, or [`chrom_peaks`](@ref) for a
trace), a `Vector{<:AbstractPeak}`, a bare m/z (or abscissa) vector, or
`mz => "text"` pairs. Useful keywords:

| Keyword | Effect |
|---------|--------|
| `nlabels` | maximum number of labels (default 10) |
| `tol` / `ppm` | apex-match tolerance for the supplied m/z |
| `declutter` | `:suppress` (default), `:none`, or a minimum gap in m/z |
| `show_mz` | append the matched m/z to a text label |
| `show_charge` | append a charge glyph (e.g. `2+`) for multiply-charged ions |
| `annot_fmt` | `printf` format for the m/z (default `"%.2f"`) |
| `color_by` | `:none` (default) or `:group` (colour annotations by group) |

Under an interactive backend (`plotly()`), each labelled peak also carries a
hover tooltip with its m/z, charge, and label. The keyword is named `annot` (not
`annotate`) because `Plots` reserves `annotate` as an alias of its own
`annotations` attribute.

## Peak-window overlays

Where `annot` *labels* peaks, a [`MassJ.Peak`](@ref) or [`MassJ.TargetPeak`](@ref)
can be **overlaid** on a spectrum with `plot!`, shading its integration window(s)
as vertical bands — so you can see exactly where a peak integrates (for a cluster
or formula `TargetPeak`, every isotopologue sub-window is drawn):

```julia
plot(avg)                                                     # the spectrum
plot!(TargetPeak("Nd(NO3)4", "Precursor"; charge = -1, tol = 0.2))
plot!(peaks)                                                  # a whole list at once
plot!(tp; seriescolor = :red, seriesalpha = 0.3, label = "window")
```

`plot!(peaks)` cycles a colour per window and labels each by `peak.label`. Combine
the two layers to both shade the window and label the apex:

```julia
plot(avg, annot = [tp]); plot!(tp)
```

For a [`MassJ.YieldCurve`](@ref), the x-axis is taken from `yc.x` and labelled with `yc.xlabel`, and one line is drawn per peak using `yc.labels` for the legend. When `yc.yields_err` contains any finite values, 1-σ ribbons are drawn around each line with `fillalpha = 0.15`. Either default can be overridden by passing the same keyword to `plot`:

```julia
plot(yc)                                   # ribbon if errors are available
plot(yc; ribbon = nothing)                 # turn ribbon off
plot(yc; fillalpha = 0.30)                 # darker ribbon
```

To drop one or more peaks from a plot — typically a dominant precursor that
swamps the axes — use [`drop_peaks`](@ref):

```julia
plot(drop_peaks(yc, "precursor"))
plot(drop_peaks(yc, ["precursor", "solvent"]))
```

## Covariance / association maps

A square association map (from [`covariance_matrix`](@ref),
[`partial_correlation`](@ref), or [`cmi_matrix`](@ref)) is drawn as a diverging
heatmap centred on zero, so positive (same-precursor) and negative
(different-precursor) associations read symmetrically:

```julia
using MassJ.plots, Plots
A = covariance_matrix(am.matrix)
covmap(cut_autocorrelation(A), am.mz)       # heatmap; suppress the diagonal first
```

[`MassJ.plots.covmap_marginal`](@ref) draws the classic covariance-mapping figure —
the map as a central heatmap with the 1-D marginal spectrum above (shared m/z on the
x-axis) and to the left (shared m/z on the y-axis), so a cross-peak lines up with the
two ions that produced it:

```julia
covmap_marginal(A, am.mz, vec(sum(am.matrix, dims = 1)))
```

Both work with any Plots backend — `gr()` for print figures, `plotly()` for
interactive notebooks. See [Chimeric tandem mass spectra](@ref) for the analysis that
produces the map and its picked pairs.
