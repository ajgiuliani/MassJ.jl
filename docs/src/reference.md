```@meta
CurrentModule = MassJ
DocTestSetup  = quote
    using LightXML
end
```

This page lists the public API of the `MassJ.jl` package. Internal helpers
that implement these functions are intentionally not exposed here; their
docstrings live in the source and are reachable via the REPL `?name` help.

```@contents
Pages = ["reference.md"]
```

# Main module
```@docs
MassJ
```

## Types
--------

### Data types
```@docs
MassJ.MScontainer
MassJ.MSscan
MassJ.MSscans
MassJ.Chromatogram
MassJ.Mobilogram
MassJ.Ionogram
MassJ.Isotope
MassJ.AbstractPeak
MassJ.Peak
MassJ.TargetPeak
MassJ.YieldCurve
```

### Method types
```@docs
MassJ.MethodType
MassJ.BasePeak
MassJ.TIC
MassJ.∆MZ
MassJ.MZ
MassJ.SG
MassJ.TBPD
MassJ.SNRA
MassJ.TopHat
MassJ.LOESS
MassJ.IPSA
MassJ.UniDec
MassJ.Charges
MassJ.Masses
```

### Filters
```@docs
MassJ.FilterType
MassJ.RT
MassJ.IC
MassJ.Level
MassJ.Scan
MassJ.Polarity
MassJ.Activation_Method
MassJ.Activation_Energy
MassJ.Precursor
MassJ.DriftTime
MassJ.CompensationVoltage
```

## I/O
------
Top-level entry points for reading and summarising spectrum files. The
appropriate format reader is selected automatically from the file extension
(mzXML, mzML, MGF, MSP, imzML, TXT).

```@docs
MassJ.info(filename::String; verbose::Bool = false)
MassJ.load(filename::String)
MassJ.retention_time(filename::String)
MassJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
MassJ.average(filename::String, arguments::FilterType...; stats::Bool=true)
```

## Filtering on loaded scans
The same operations are available directly on a `Vector{MSscan}` produced by
[`load`](@ref). Filters are composed into a single short-circuiting predicate
and applied in one pass — see [Composed predicates](@ref).

```@docs
MassJ.average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
MassJ.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
MassJ.retention_time(scans::Vector{MSscan})
```

## Extracting subsets
```@docs
MassJ.extract(filename::String, arguments::FilterType...)
MassJ.extract(scans::Vector{MSscan}, arguments::FilterType...)
```

## Composed predicates
Internal predicate composition used by `extract`/`chromatogram`/`average` on a
`Vector{MSscan}`. Exposed here for users extending the package with new
[`MassJ.FilterType`](@ref) subtypes — implement `to_predicate(f::MyFilter)`
and the new filter is automatically composable with all existing ones.

```@docs
MassJ.to_predicate
MassJ.compose_predicates
```

## Processing
-------------

### Mass spectrum
```@docs
MassJ.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
MassJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
MassJ.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
MassJ.centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100))
MassJ.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
MassJ.baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
```

## Deconvolution
----------------
```@docs
MassJ.deconv
```

## Simulations
--------------
```@docs
MassJ.formula
MassJ.masses
MassJ.isotopic_distribution
MassJ.simulate
```


## Energy-resolved yields
-------------------------
```@docs
MassJ.yields
MassJ.integrate_window
MassJ.normalize_tic
MassJ.normalize_flux
MassJ.drop_peaks
MassJ.read_peaklist
MassJ.write_csv
```


# Plots
-------
```@autodocs
Modules = [MassJ.plots]
```


# Utilities
-----------

## Base overloaded operators
Arithmetic is overloaded on [`MassJ.MScontainer`](@ref) so that scans can be
combined directly — interpolation aligns the m/z axes automatically.

```@docs
+(a::MScontainer, b::MScontainer)
-(a::MScontainer, b::MScontainer)
/(a::MSscan, N::Real)
/(a::MSscans, N::Real)
*(a::MSscan, N::Real)
*(a::MSscans, N::Real)
*(N::Real, a::MScontainer)
*(a::MScontainer, b::MScontainer)
```

## Other utilities
```@docs
MassJ.avg(a::MScontainer, b::MScontainer)
MassJ.num2pnt(x::Vector{Float64}, val::Real)
```
