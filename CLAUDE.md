# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MassJ.jl is a Julia package for loading, processing, and analyzing mass spectrometry data. It supports mzXML and TXT file formats and provides signal processing, isotopic distribution calculations, charge/mass deconvolution, and plotting capabilities.

## Common Commands

```bash
# Run tests (from repo root)
cd test && julia runtests.jl

# Run tests via Pkg
julia -e 'using Pkg; Pkg.test("MassJ")'

# Build documentation
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl

# Install package locally for development
julia -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd()))'
```

## Architecture

### Module structure (`src/MassJ.jl`)

The main module includes files in a specific order — **ordering matters** because later files depend on types/functions from earlier ones:

1. **types.jl** — Abstract type hierarchy and core structs
2. **Io.jl** — Public API for file I/O (`load`, `info`, `chromatogram`, `average`, `extract`, `retention_time`)
3. **msscans.jl** — Internal filtering/dispatch by `FilterType` subtypes
4. **mzxml.jl** — mzXML parser (LightXML + zlib decompression)
5. **process.jl** — Signal processing: `smooth`, `centroid`, `baseline_correction`
6. **extract.jl** — Data extraction by filter criteria
7. **utilities.jl** — Overloaded `+`, `-`, `*`, `/` operators on `MScontainer` types
8. **isotopes.jl** + isotopes-data.jl — Formula parsing, mass calculation, isotopic distributions
9. **deconvolution.jl** — Charge/mass deconvolution
10. **txt.jl** — Plain text file loader
11. **plots.jl** — `MassJ.plots` submodule with Plots.jl recipes

### Type system (dispatch-driven design)

The package relies heavily on Julia's multiple dispatch. There are three abstract type hierarchies:

- **`MScontainer`** — Data containers: `MSscan` (single scan), `MSscans` (averaged scans), `Chromatogram`
- **`MethodType`** — Algorithm selectors passed to processing functions (e.g., `SG()` for smoothing, `TopHat()` for baseline, `TBPD()` for centroiding, `TIC()`/`BasePeak()` for chromatograms)
- **`FilterType`** — Data selectors passed to `chromatogram`/`average`/`extract` (e.g., `RT([start,end])`, `Level(n)`, `Precursor(mz)`)

Processing functions like `smooth(scan; method=SG(...))` dispatch on the `MethodType` to select the algorithm.

### Operators

`Base.+`, `Base.-`, `Base.*`, `Base./` are overloaded in `utilities.jl` for `MScontainer` types. Arithmetic between two scans performs interpolation to align m/z axes before combining.

## Test Data

Test files in `test/`: `test.mzXML` (6 scans, mixed MS levels), `test64.mzXML` (64-bit encoded), `bad1-3.mzXML` (malformed files for error handling tests).
