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
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ajgiuliani.github.io/MassJ.jl/dev)

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
Documentation is available [here](https://ajgiuliani.github.io/MassJ.jl/dev).


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

See the [documentation](https://ajgiuliani.github.io/MassJ.jl/dev) for additional information.

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

## Changelog

### v1.0.0 (2026-04-09)

**Major release — Package renamed from MSj.jl to MassJ.jl**

#### Breaking Changes
- Renamed package from `MSj` to `MassJ` (all imports must be updated)

#### New Features
- Added mzML file format reader
- Added MGF (Mascot Generic Format) file format reader
- Added MSP file format reader
- Added imzML (imaging mass spectrometry) file format reader
- Added ion mobility data types support
- Rewrote charge/mass deconvolution algorithm
- Optimized isotopic distribution calculations

#### Improvements
- Switched CI from Travis CI to GitHub Actions with coverage
- Switched docs deployment to GitHub Actions with modern Documenter.jl
- Bumped minimum Julia version to 1.10 (LTS), tested on 1.12
- Updated dependencies: fixed Interpolations `Line` import, updated LsqFit compat
- Added AppVeyor CI configuration
- Improved float precision handling in tests for nightly Julia
- Fixed docstring typos in mzxml.jl
- Updated documentation for deconvolution, isotopes, and public API
- Added tests for isotopes, deconvolution, and interpolation import
- Minor cleanup in `utilities.jl` (convolve, same-size array handling)
- Minor cleanup in `process.jl`

#### License
- Changed license to CeCILL-2.1, then to GPL-3.0

---

### v0.3.1 (2019-12-12)

#### New Features
- Added zlib compression support for mzXML files (via Libz dependency)

#### Bug Fixes
- Fixed bug in mzxml.jl `ActivationMethod` parsing (added conditional check)
- Fixed bug in TopHat baseline filtering
- Fixed `Int64` to `Int` for cross-platform compatibility

#### Improvements
- Renamed internal module from `msJ` to `MSJ`/`MSj` (package identity consolidation)
- Added CompatHelper for automated dependency version bumping
- Added AppVeyor CI configuration
- Updated isotopes: fixed formula parsing, changed trim conditions so that `p >= target`
- Documentation updates and notebook fixes
- Updated CI to Julia 1.3

---

### v0.3 (2019-11-21)

#### New Features
- **Isotopic distribution calculations** — new `isotopes.jl` module with:
  - Chemical formula parsing
  - Molecular mass calculation
  - Isotopic distribution simulation

#### Improvements
- Removed DSP.jl dependency — replaced with internal convolution function in `utilities.jl`
- Updated dependency versions
- Added notebook tutorial for isotopic distributions
- Documentation updates and typo fixes

---

### v0.2.1 (2019-10-23)

#### Improvements
- Fixed documentation URLs in README
- Updated and fixed Jupyter notebook examples (examples 2 and 3)
- Added Documenter.jl compat entry
- Tutorial and documentation typo fixes

---

### v0.2 (2019-10-21)

#### New Features
- **Baseline correction** — multiple algorithms:
  - `TopHat()` — morphological top-hat filtering (via ImageMorphology)
  - `LOESS()` — locally estimated scatterplot smoothing
  - `IPSA()` — iterative polynomial smoothing (default method)
- **Extract function** — new `extract.jl` for extracting subsets of scan data by index
- **SNR peak picking** — added signal-to-noise ratio method to `centroid`
- **Morphological filtering** — added morphological operations in `utilities.jl`
- Made `average` function public (renamed from `msfilter`)

#### Bug Fixes
- Fixed bug in pseudo-Voigt profile calculation
- Fixed centroid and baseline_correction failing tests
- Fixed extract bug when index = 1
- Fixed `basePeakIntensity` calculation after smoothing

#### Improvements
- Added Lorentz and Voigt peak profiles to `process.jl`
- Simplified Savitzky-Golay function calls
- Updated centroid default parameters
- Added `MSscan`/`MSscans` type descriptions in documentation
- Added smooth and centroid processing for `Vector{MSscan}`
- Updated CI to Julia 1.3 (removed 1.2)
- Added tutorials with Jupyter notebooks
- Extensive documentation updates

---

### v0.1.1 (2019-05-30)

**Initial registered release**

#### Features
- mzXML file loading and parsing (LightXML-based)
- Single scan and averaged scan data containers (`MSscan`, `MSscans`)
- Chromatogram extraction (`TIC`, `BasePeak`)
- Signal processing: `smooth` (Savitzky-Golay), `centroid` (TBPD)
- Arithmetic operators (`+`, `-`, `*`, `/`) on mass spec containers with m/z interpolation
- Charge/mass deconvolution
- Plot recipes via Plots.jl/RecipesBase
- Travis CI integration
- Basic documentation with Documenter.jl
