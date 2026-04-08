```@meta
CurrentModule = MassJ
DocTestSetup  = quote
    using LightXML
end
```

This page lists all the documented elements of the `MassJ.jl` package covering all modules and submodules.

```@contents
Pages = ["reference.md"]
```

# Main module
```@docs
MassJ
```

## Types
--------
Submodule with types and structures used to stored the data and dispatch to the right methods.

### Data types
```@docs
MassJ.MScontainer
MassJ.MSscan
MassJ.MSscans
MassJ.Chromatogram
MassJ.Mobilogram
MassJ.Ionogram
MassJ.Isotope
```

### Methods Types
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

Module for importing and exporting data. Dispatch to specific methods according to the file extension

```@docs
MassJ.info(filename::String; verbose::Bool = false)
MassJ.load(filename::String)
MassJ.retention_time(filename::String)
MassJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
MassJ.average(filename::String, arguments::FilterType...; stats::Bool=true)
```

### mzXML
---------
Interface to the mzXML file format.

```@docs
MassJ.info_mzxml
MassJ.load_mzxml_all
MassJ.load_mzxml
MassJ.load_mzxml_spectrum
MassJ.retention_time(msRun::XMLElement)
MassJ.filter(msRun::XMLElement, argument::Level{<:Int})
MassJ.filter(msRun::XMLElement, argument::Level{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::Scan{<:Int})
MassJ.filter(msRun::XMLElement, argument::Scan{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::Polarity{<:String})
MassJ.filter(msRun::XMLElement, argument::Polarity{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::RT{<:Real})
MassJ.filter(msRun::XMLElement, argument::RT{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
MassJ.filter(msRun::XMLElement, argument::IC{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::Precursor{<:Real})
MassJ.filter(msRun::XMLElement, argument::Precursor{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::Activation_Energy{<:Real})
MassJ.filter(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
MassJ.filter(msRun::XMLElement, argument::Activation_Method{<:String})
MassJ.filter(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
MassJ.extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
MassJ.composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
```

### mzML
---------
Interface to the mzML file format (PSI standard). Uses the same binary decoding pipeline as mzXML (Codecs, Libz) but with little-endian byte order and separate m/z and intensity arrays.

```@docs
MassJ.info_mzml
MassJ.load_mzml_all
MassJ.load_mzml
MassJ.load_mzml_spectrum
MassJ.retention_time_mzml
MassJ.find_mzml_root
MassJ.get_cv_param
MassJ.get_cv_value
MassJ.has_cv_param
```

### MGF
--------
Interface to the MGF (Mascot Generic Format) file format. Text-based format for centroided MS/MS peak lists.

```@docs
MassJ.load_mgf_all
MassJ.build_mgf_scan
MassJ.info_mgf
```

### MSP
--------
Interface to the MSP (NIST Mass Spectral Library) file format. Text-based format used by NIST, MoNA, MassBank, and GNPS for spectral libraries.

```@docs
MassJ.load_msp_all
MassJ.build_msp_scan
MassJ.info_msp
```

### imzML
----------
Interface to the imzML file format for imaging mass spectrometry. Consists of an XML metadata file (`.imzML`) and a companion binary data file (`.ibd`). Supports both continuous and processed storage modes, and stores spatial coordinates in the metadata dictionary.

```@docs
MassJ.load_imzml_all
MassJ.load_imzml_spectrum
MassJ.info_imzml
```


## Filtering
```@docs
MassJ.average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
MassJ.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
MassJ.retention_time(scans::Vector{MSscan})
MassJ.filter(scans::Vector{MSscan}, argument::Scan{<:Int})
MassJ.filter(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::Level{<:Int})
MassJ.filter(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::Precursor{<:Real})
MassJ.filter(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
MassJ.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::Activation_Method{<:String})
MassJ.filter(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::Polarity{<:String})
MassJ.filter(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::RT{<:Real}) 
MassJ.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
MassJ.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
MassJ.filter(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
MassJ.extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
MassJ.composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
```

## Extracting subsets
```@docs
MassJ.extract(filename::String, arguments::FilterType...)
MassJ.extract(scans::Vector{MSscan}, arguments::FilterType...)
MassJ.build_subset(filename::String, indices::Vector{Int})
MassJ.build_subset(scans::Vector{MSscan}, indices::Vector{Int})
```

## Process
----------

### Mass spectrum
-----------------
```@docs
MassJ.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
MassJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
MassJ.savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
MassJ.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
MassJ.centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100))
MassJ.snra(scan::MScontainer, thres::Real, region::Int)
MassJ.tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)
MassJ.gauss(x::Float64, p::Vector{Float64})
MassJ.lorentz(x::Float64, p::Vector{Float64})
MassJ.voigt(x::Float64, p::Vector{Float64})
MassJ.tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)
MassJ.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
MassJ.baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
MassJ.tophat_filter(scan::MScontainer, region::Int )
MassJ.tophat_filter(scans::Vector{MSscan}, region::Int )
MassJ.loess(scans::Vector{MSscan}, iter::Int )
MassJ.loess(scan::MScontainer, iter::Int )
MassJ.ipsa(scan::MScontainer, width::Real, maxiter::Int)
MassJ.ipsa(scans::Vector{MSscan}, width::Real, maxiter::Int)
```


### Chromatogram
----------------
No functions yet. To be added.


## Deconvolution
----------------
```@docs
MassJ.deconv
```

## Simulations
-------------
```@docs
MassJ.formula
MassJ.masses
MassJ.isotopic_distribution
MassJ.simulate
```


# Plots
--------
```@autodocs
Modules = [MassJ.plots]

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
MassJ.avg(a::MScontainer, b::MScontainer)
MassJ.add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
MassJ.num2pnt(x::Vector{Float64}, val::Real)
MassJ.savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
MassJ.extremefilt(input::AbstractArray, minmax::Function, region::Int)
MassJ.morpholaplace(input::AbstractArray, region::Int)
MassJ.morphogradient(input::AbstractArray, region::Int)
MassJ.tophat(input::AbstractArray, region::Int)
MassJ.bottomhat(input::AbstractArray, region::Int) 
MassJ.opening(input::AbstractArray, region::Int)
MassJ.closing(input::AbstractArray, region::Int)
MassJ.erosion(input::AbstractArray, region::Int)
MassJ.dilatation(input::AbstractArray, region::Int)
MassJ.convolve(a::AbstractArray, b::AbstractArray)
```
