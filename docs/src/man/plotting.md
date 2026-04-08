# Plotting
Plotting facilities are available as a submodule to the `MassJ` package.  The [`MassJ.plots`](@ref) module relies on the [RecipesBase package](https://github.com/JuliaPlots/RecipesBase.jl), which allows writing recipes to plot users' data types. Hence, recipes have been created for `MSscan`, `Msscans` and `Chromatogram`:

```julia
plot(scans[1], method = :relative))

```

By default plotting is made in relative intensities, which may be changed by setting method to `:absolute`.
