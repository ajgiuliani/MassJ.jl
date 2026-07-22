```@meta
CurrentModule = MassJ
DocTestSetup  = quote
    using LightXML
end
```

# References

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
MassJ.MSscans
MassJ.IonCurrent
MassJ.maxic
MassJ.MSrun
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
MassJ.CWT
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
MassJ.ChargeState
MassJ.SpectrumType
MassJ.MobilityType
MassJ.MetaData
MassJ.InstrumentConfig
```

## I/O
------
Top-level entry points for reading and summarising spectrum files. The
appropriate format reader is selected automatically from the file extension
(mzXML, mzML, MGF, MSP, imzML, TXT).

```@docs
MassJ.info(filename::String; verbose::Bool = false)
MassJ.load(filename::String)
MassJ.load(paths::AbstractVector{<:AbstractString})
MassJ.load_usi
MassJ.retention_time(filename::String)
MassJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
MassJ.mobilogram
MassJ.ionogram
MassJ.average(filename::String, arguments::FilterType...; stats::Bool=true)
MassJ.average(paths::AbstractVector{<:AbstractString}, arguments::FilterType...; stats::Bool=true, recursive::Bool=false)
MassJ.save
MassJ.save_mzml
MassJ.save_mzxml
MassJ.save_mgf
MassJ.save_msp
MassJ.save_txt
MassJ.save_imzml
```

## Filtering on loaded scans
The same operations are available directly on the `MSrun` (or any
`AbstractVector{MSscans}`) produced by [`load`](@ref). Filters are composed into
a single short-circuiting predicate and applied in one pass — see
[Composed predicates](@ref).

```@docs
MassJ.average(scans::AbstractVector{MSscans}, arguments::FilterType...; stats::Bool=true)
MassJ.chromatogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
MassJ.retention_time(scans::AbstractVector{MSscans})
```

## Extracting subsets
```@docs
MassJ.extract(filename::String, arguments::FilterType...)
MassJ.extract(scans::AbstractVector{MSscans}, arguments::FilterType...)
```

## Composed predicates
Internal predicate composition used by `extract`/`chromatogram`/`average` on any
`AbstractVector{MSscans}`. Exposed here for users extending the package with new
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
MassJ.smooth(scans::AbstractVector{MSscans}; method::MethodType=SG(5, 9, 0))
MassJ.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
MassJ.centroid(scans::AbstractVector{MSscans}; method::MethodType=SNRA(1., 100))
MassJ.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
MassJ.baseline_correction(scans::AbstractVector{MSscans}; method::MethodType=TopHat(100) )
```

### Ion-current trace (chromatogram / mobilogram / ionogram)
```@docs
MassJ.smooth(ic::IonCurrent; method::MethodType=SG(5, 9, 0))
MassJ.baseline_correction(ic::IonCurrent; method::MethodType=TopHat(100))
MassJ.ChromPeak
MassJ.chrom_peaks
MassJ.chrom_peak
```

### 2-D LC-MS feature detection
Detect m/z × retention-time features in a single run: build mass traces, then run
chromatographic peak detection on each.
```@docs
MassJ.Feature
MassJ.detect_features
MassJ.feature_table
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
MassJ.isotope_table
MassJ.simulate
```

## Adducts and neutral-mass conversion
```@docs
MassJ.Adduct
MassJ.Adducts
MassJ.adduct_mz
MassJ.neutral_mass
```

## Mass calibration
```@docs
MassJ.Calibration
MassJ.calibrate
```

## Formula assignment from m/z
```@docs
MassJ.FormulaCandidate
MassJ.assign_formula
MassJ.score_isotope_pattern
```

## Energy-resolved yields
-------------------------
```@docs
MassJ.yields
MassJ.integrate_window
MassJ.normalize_tic
MassJ.normalize_external
MassJ.normalize_flux
MassJ.combine_yields
MassJ.shift_x
MassJ.scale_yields
MassJ.recalibrate_x
MassJ.trim_yields
MassJ.restrict_x
MassJ.drop_peaks
MassJ.read_peaklist
MassJ.write_csv
```


## Chimeric-spectra statistical decomposition
----------------------------------------------
Data-driven resolution of multiplexed (chimeric) tandem mass spectra. The score
matrices feed hierarchical clustering to group ions by precursor.
[`MassJ.cmi_matrix`](@ref) requires `EntropyInvariant` and
[`MassJ.cluster_ions`](@ref) requires `Clustering`; both are loaded on demand as
package extensions.

Each association can carry a per-edge significance — a leave-one-out
[`MassJ.jackknife_significance`](@ref) z-score for the covariance/correlation family,
or a [`MassJ.permutation_significance`](@ref) p-value for conditional mutual
information — which [`MassJ.fdr_adjust`](@ref) turns into a Benjamini-Hochberg gate.
[`MassJ.correspondence`](@ref) chains the whole pipeline (abundance matrix → score →
significance → FDR gate → clustering → per-precursor spectra) into one hands-off call.

```@docs
MassJ.abundance_matrix
MassJ.partial_correlation
MassJ.covariance_matrix
MassJ.cmi_matrix
MassJ.jackknife_significance
MassJ.permutation_significance
MassJ.fdr_adjust
MassJ.cluster_ions
MassJ.cluster_spectra
MassJ.correspondence
```


## Covariance maps and fragment-pair features
---------------------------------------------
The pairwise route, complementary to the clustering above: pick correlated
fragment-ion *pairs* off a two-dimensional association map, integrate each pair's
covariance volume, and rate it by a leave-one-out jackknife signal-to-noise. Plot the
map with [`MassJ.plots.covmap`](@ref) or the framed figure
[`MassJ.plots.covmap_marginal`](@ref).

```@docs
MassJ.cut_autocorrelation
MassJ.MapFeature
MassJ.map_features
MassJ.covariance_features
```


## Ecosystem interoperability
----------------------------
Bridges between MassJ's typed containers and the wider Julia data / machine-learning
ecosystem: featurise a spectrum series to a matrix, map results back, and rebuild
spectra. See [Interoperability](@ref).

```@docs
MassJ.featurize
MassJ.select_spectra
MassJ.from_matrix
MassJ.spectra_table
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
/(a::MSscans, N::Real)
*(a::MSscans, N::Real)
*(N::Real, a::MScontainer)
*(a::MScontainer, b::MScontainer)
```

## Uncertainty accessors
Retrieve the per-m/z uncertainty carried by an averaged spectrum — see
[Uncertainties and errors](@ref).

```@docs
MassJ.nscans
MassJ.stdev
MassJ.sem
MassJ.measurements
MassJ.withunits
```

## Other utilities
```@docs
MassJ.avg(a::MScontainer, b::MScontainer)
MassJ.num2pnt(x::Vector{Float64}, val::Real)
```
