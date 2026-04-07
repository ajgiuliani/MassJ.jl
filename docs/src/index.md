Home
====

`MSj.jl` is a package for loading, processing and plotting mass spectrometry data. It provides a range of functionalities such as:
* Getting information on the file
* Load a file
* Averaging mass spectra based on various criteria that may be combined
* Chromatogram
* Processing the data
   - smoothing
   - baseline correction
   - peak-picking
* Calculation of isotopic distribution

The [tutorial page](tutos/tutorials.md) shows examples how to use this package and provides a general guide to it. The [manual page](man/introduction.md) explains the structure of the package and the [reference page](reference.md) gives a full documentation for each type and function.


## Installation
---------------
There are two ways of installing the package.
```julia
julia> using Pkg ;
julia> Pkg.add(PackageSpec(url="https://github.com/ajgiuliani/MSj.jl"))
```

or

```julia
julia> ]
(v1.3) pkg>  add https://github.com/ajgiuliani/MSj.jl
```


## Supported file formats
-----------
  - mzXML
  - mzML (PSI standard)
  - MGF (Mascot Generic Format)
  - MSP (NIST Mass Spectral Library)
  - imzML (Imaging Mass Spectrometry)
  - TXT (plain text peak lists)
