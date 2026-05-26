<p align="center">
  <h1 align="center"> MassJ.jl  </h1>
</p>
<p align="center">
  <img align="center" src="docs/src/assets/logo.png" width="400" height="200" />
</p>
<p align="center">
  <normal> A mass spectrometry package for Julia </normal>
</p>



##
[![CI](https://github.com/ajgiuliani/MassJ.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ajgiuliani/MassJ.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/ajgiuliani/MassJ.jl/badge.svg?branch=master)](https://coveralls.io/github/ajgiuliani/MassJ.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ajgiuliani.github.io/MassJ.jl/stable)

## Installation
This package is registered. It can be installed either with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```julia
pkg> add MassJ
```
or using the package API:

```julia
using Pkg
Pkg.add("MassJ")
```

## Documentation
Documentation is available [here](https://ajgiuliani.github.io/MassJ.jl/stable).


## Usage
MassJ is a package for loading, processing and plotting mass spectrometry data. It provides the following functionalities:

    Getting information on the file
    Load a file
    Averaging mass spectra based on various criteria that may be combined
    Chromatogram and extracted chromatograms
    Processing the data
        smoothing
        baseline correction
        peak-picking
    Calculation of isotopic distribution
	Charge states deconvolution

To get information on a file:
```julia
info("path/to/file")
```

Mass spectra can be loaded by:
```julia
data = load("path/to/file")
```

And averaged as follow:
```julia
ms1 = average(data, MassJ.Level(1))                   # full MS scans
ms2 = average(data, MassJ.Level(2))                   # MS2 spectra
ms3 = average(data, MassJ.Activation_Method("CID"))   # CID spectra
```

See the [documentation](https://ajgiuliani.github.io/MassJ.jl/stable) for additional information.

## Supported file format
* mzxml
* mzML
* MGF
* MSP
* imzML
* txt

## Other Julia packages
* [MzXML.jl](https://github.com/timholy/MzXML.jl): Load mass spectrometry mzXML files.
* [MassSpec.jl](https://gitlab.com/odurif/MassSpec.jl): Mass spectometry utilities for Julia

## References
The isotopic distribution and charge deconvolution routines are independent reimplementations of published algorithms. If you use them, please cite the original works:

* **IsoSpec** (isotopic distributions): Łącki, M. K.; Startek, M.; Valkenborg, D.; Gambin, A. *IsoSpec: Hyperfast Fine Structure Calculator.* Anal. Chem. **2017**, 89, 3272. [Project](https://github.com/MatteoLacki/IsoSpec)
* **UniDec** (charge/mass deconvolution): Marty, M. T.; Baldwin, A. J.; Marklund, E. G.; Hochberg, G. K. A.; Benesch, J. L. P.; Robinson, C. V. *Bayesian Deconvolution of Mass and Ion Mobility Spectra: From Binary Interactions to Polydisperse Ensembles.* Anal. Chem. **2015**, 87, 4370. DOI: [10.1021/acs.analchem.5b00140](https://doi.org/10.1021/acs.analchem.5b00140). [Project](https://github.com/michaelmarty/UniDec)

