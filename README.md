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
[![codecov](https://codecov.io/gh/ajgiuliani/MassJ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ajgiuliani/MassJ.jl)
[![Coverage Status](https://coveralls.io/repos/github/ajgiuliani/MassJ.jl/badge.svg?branch=master)](https://coveralls.io/github/ajgiuliani/MassJ.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ajgiuliani.github.io/MassJ.jl/stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## What MassJ does

MassJ.jl is a general-purpose Julia package for loading, processing, and
analysing mass spectromet8ry data. It is designed to slot into the rest of the
Julia ecosystem (Plots.jl, Statistics, DataFrames) rather than to replicate
the feature surface of pipelines like OpenMS or pyteomics.

- **Read & write open formats** — mzML (incl. indexed-mzML), mzXML, MGF, MSP,
  imzML, plain TXT. mzML output is PSI 1.1.0 schema-validated and accepted by
  MaxQuant and other downstream tools.
- **Signal processing** — Savitzky-Golay smoothing, centroiding (SNRA,
  template-based peak detection, continuous wavelet transform), baseline
  correction (TopHat, LOESS, IPSA).
- **Isotopes & deconvolution** — chemical-formula parsing, exact/average/
  nominal masses, isotopic distribution simulation, UniDec-style charge-state
  deconvolution.
- **Energy-resolved yields** — purpose-built pipeline for UVPD / action
  spectroscopy / CID breakdown: locate peaks, integrate per scan, propagate
  proper per-scan errors, normalize by TIC and/or photon flux, plot with
  uncertainty ribbons.
- **Native plots** — Plots.jl recipes for spectra, chromatograms, and yield
  curves out of the box.

## Installation

```julia
julia> ]
pkg> add MassJ
```

## Quick start

```julia
using MassJ, Plots

scans = load("input.mzML")               # → MSrun (a Vector{MSscan} + file metadata)
avg   = average(scans, Level(1))         # average all MS1 spectra
clean = baseline_correction(smooth(avg)) # SG smooth + TopHat baseline
plot(clean)                              # one-line plot

# Save back (indexed mzML by default; MaxQuant-compatible)
save(clean, "processed.mzML")
```

A more substantial workflow — energy-resolved yields for action-spectroscopy
data:

```julia
peaks = [Peak(110.5, "fragment_a"; tol = 0.5),
         Peak(500.05, "precursor"; ppm = 5.0)]

yc = yields("data/UVPD/", peaks; x0 = 3.0, step = 0.1,
            xlabel = "photon energy (eV)")
yc = yc |> normalize_tic |> y -> normalize_flux(y, "flux.txt")

plot(yc)                                 # ribbons show per-scan SEM
write_csv(yc, "yields.csv")
```

## Documentation

Full documentation is at
[ajgiuliani.github.io/MassJ.jl](https://ajgiuliani.github.io/MassJ.jl/stable).
The manual covers each format, the processing pipeline, the yields workflow,
and the reference API.

## Supported file formats

| Format | Read | Write          |
|--------|:----:|:--------------:|
| mzML   | ✅   | ✅ (indexed)   |
| mzXML  | ✅   | ✅              |
| MGF    | ✅   | ✅              |
| MSP    | ✅   | ✅              |
| imzML  | ✅   | ✅              |
| TXT    | ✅   | ✅              |

## Citing

If you use MassJ in published work, please cite it via the [CITATION.cff](CITATION.cff) metadata at the repository root.

The deconvolution, isotopic-distribution, and wavelet peak-detection routines
are independent re-implementations of published algorithms — please also cite
the original work:

* **IsoSpec** (isotopic distributions): Łącki, M. K.; Startek, M.;
  Valkenborg, D.; Gambin, A. *IsoSpec: Hyperfast Fine Structure Calculator.*
  Anal. Chem. **2017**, 89, 3272.
  [Project](https://github.com/MatteoLacki/IsoSpec)
* **CWT peak detection** (continuous wavelet transform centroiding): Du, P.;
  Kibbe, W. A.; Lin, S. M. *Improved peak detection in mass spectrum by
  incorporating continuous wavelet transform-based pattern matching.*
  Bioinformatics **2006**, 22, 2059.
  DOI: [10.1093/bioinformatics/btl355](https://doi.org/10.1093/bioinformatics/btl355).
* **UniDec** (charge/mass deconvolution): Marty, M. T.; Baldwin, A. J.;
  Marklund, E. G.; Hochberg, G. K. A.; Benesch, J. L. P.; Robinson, C. V.
  *Bayesian Deconvolution of Mass and Ion Mobility Spectra: From Binary
  Interactions to Polydisperse Ensembles.* Anal. Chem. **2015**, 87, 4370.
  DOI: [10.1021/acs.analchem.5b00140](https://doi.org/10.1021/acs.analchem.5b00140).
  [Project](https://github.com/michaelmarty/UniDec)

## Related Julia packages

* [MzXML.jl](https://github.com/timholy/MzXML.jl) — older mzXML loader
* [julia_mzML_imzML](https://codeberg.org/LabABI/julia_mzML_imzML) — mzML files loader
* [MassSpec.jl](https://gitlab.com/odurif/MassSpec.jl) — assorted MS utilities

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md)
for the workflow.

## License

GPL-3.0-or-later, see [LICENSE](LICENSE).
