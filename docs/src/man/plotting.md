# Plotting
Plotting facilities are available as a submodule to the `MassJ` package.  The [`MassJ.plots`](@ref) module relies on the [RecipesBase package](https://github.com/JuliaPlots/RecipesBase.jl), which allows writing recipes to plot users' data types. Recipes are provided for `MSscan`, `MSscans`, `Chromatogram`, and `YieldCurve`:

```julia
plot(scans[1], method = :relative)
plot(yc)                                   # YieldCurve: one line per peak
```

For mass spectra and chromatograms, plotting is made in relative intensities by default; this can be changed by setting `method = :absolute`. For a [`MassJ.YieldCurve`](@ref), the x-axis is taken from `yc.x` and labelled with `yc.xlabel`, and one line is drawn per peak using `yc.labels` for the legend.
