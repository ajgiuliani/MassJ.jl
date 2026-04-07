```@meta
CurrentModule = MSj
DocTestSetup  = quote
    using LightXML
end
```

This page lists all the documented elements of the `MSj.jl` package covering all modules and submodules.

```@contents
Pages = ["reference.md"]
```

# Main module
```@docs
MSj
```

## Types
--------
Submodule with types and structures used to stored the data and dispatch to the right methods.

### Data types
```@docs
MSj.MScontainer
MSj.MSscan
MSj.MSscans
MSj.Chromatogram
MSj.Mobilogram
MSj.Ionogram
MSj.Isotope
```

### Methods Types
```@docs
MSj.MethodType
MSj.BasePeak
MSj.TIC
MSj.∆MZ
MSj.MZ
MSj.SG
MSj.TBPD
MSj.SNRA
MSj.TopHat
MSj.LOESS
MSj.IPSA
MSj.UniDec
MSj.Charges
MSj.Masses
```

### Filters
```@docs
MSj.FilterType
MSj.RT
MSj.IC
MSj.Level
MSj.Scan
MSj.Polarity
MSj.Activation_Method
MSj.Activation_Energy
MSj.Precursor
MSj.DriftTime
MSj.CompensationVoltage
```

## I/O
------

Module for importing and exporting data. Dispatch to specific methods according to the file extension

```@docs
MSj.info(filename::String; verbose::Bool = false)
MSj.load(filename::String)
MSj.retention_time(filename::String)
MSj.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
MSj.average(filename::String, arguments::FilterType...; stats::Bool=true)
```

### mzXML
---------
Interface to the mzXML file format.

```@docs
MSj.info_mzxml
MSj.load_mzxml_all
MSj.load_mzxml
MSj.load_mzxml_spectrum
MSj.retention_time(msRun::XMLElement)
MSj.filter(msRun::XMLElement, argument::Level{<:Int})
MSj.filter(msRun::XMLElement, argument::Level{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::Scan{<:Int})
MSj.filter(msRun::XMLElement, argument::Scan{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::Polarity{<:String})
MSj.filter(msRun::XMLElement, argument::Polarity{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::RT{<:Real})
MSj.filter(msRun::XMLElement, argument::RT{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
MSj.filter(msRun::XMLElement, argument::IC{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::Precursor{<:Real})
MSj.filter(msRun::XMLElement, argument::Precursor{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::Activation_Energy{<:Real})
MSj.filter(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
MSj.filter(msRun::XMLElement, argument::Activation_Method{<:String})
MSj.filter(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
MSj.extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
MSj.composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
```

### mzML
---------
Interface to the mzML file format (PSI standard). Uses the same binary decoding pipeline as mzXML (Codecs, Libz) but with little-endian byte order and separate m/z and intensity arrays.

```@docs
MSj.info_mzml
MSj.load_mzml_all
MSj.load_mzml
MSj.load_mzml_spectrum
MSj.retention_time_mzml
MSj.find_mzml_root
MSj.get_cv_param
MSj.get_cv_value
MSj.has_cv_param
```

### MGF
--------
Interface to the MGF (Mascot Generic Format) file format. Text-based format for centroided MS/MS peak lists.

```@docs
MSj.load_mgf_all
MSj.build_mgf_scan
MSj.info_mgf
```


## Filtering
```@docs
MSj.average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
MSj.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
MSj.retention_time(scans::Vector{MSscan})
MSj.filter(scans::Vector{MSscan}, argument::Scan{<:Int})
MSj.filter(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::Level{<:Int})
MSj.filter(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::Precursor{<:Real})
MSj.filter(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
MSj.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::Activation_Method{<:String})
MSj.filter(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::Polarity{<:String})
MSj.filter(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::RT{<:Real}) 
MSj.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
MSj.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
MSj.filter(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
MSj.extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
MSj.composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
```

## Extracting subsets
```@docs
MSj.extract(filename::String, arguments::FilterType...)
MSj.extract(scans::Vector{MSscan}, arguments::FilterType...)
MSj.build_subset(filename::String, indices::Vector{Int})
MSj.build_subset(scans::Vector{MSscan}, indices::Vector{Int})
```

## Process
----------

### Mass spectrum
-----------------
```@docs
MSj.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
MSj.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
MSj.savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
MSj.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
MSj.centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100))
MSj.snra(scan::MScontainer, thres::Real, region::Int)
MSj.tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)
MSj.gauss(x::Float64, p::Vector{Float64})
MSj.lorentz(x::Float64, p::Vector{Float64})
MSj.voigt(x::Float64, p::Vector{Float64})
MSj.tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)
MSj.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
MSj.baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
MSj.tophat_filter(scan::MScontainer, region::Int )
MSj.tophat_filter(scans::Vector{MSscan}, region::Int )
MSj.loess(scans::Vector{MSscan}, iter::Int )
MSj.loess(scan::MScontainer, iter::Int )
MSj.ipsa(scan::MScontainer, width::Real, maxiter::Int)
MSj.ipsa(scans::Vector{MSscan}, width::Real, maxiter::Int)
```


### Chromatogram
----------------
No functions yet. To be added.


## Deconvolution
----------------
```@docs
MSj.deconv
```

## Simulations
-------------
```@docs
MSj.formula
MSj.masses
MSj.isotopic_distribution
MSj.simulate
```


# Plots
--------
```@autodocs
Modules = [MSj.plots]

```


# Utilities
------------

## Base overloaded
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

## Utility function
```@docs
MSj.avg(a::MScontainer, b::MScontainer)
MSj.add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
MSj.num2pnt(x::Vector{Float64}, val::Real)
MSj.savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
MSj.extremefilt(input::AbstractArray, minmax::Function, region::Int)
MSj.morpholaplace(input::AbstractArray, region::Int)
MSj.morphogradient(input::AbstractArray, region::Int)
MSj.tophat(input::AbstractArray, region::Int)
MSj.bottomhat(input::AbstractArray, region::Int) 
MSj.opening(input::AbstractArray, region::Int)
MSj.closing(input::AbstractArray, region::Int)
MSj.erosion(input::AbstractArray, region::Int)
MSj.dilatation(input::AbstractArray, region::Int)
MSj.convolve(a::AbstractArray, b::AbstractArray)
```
