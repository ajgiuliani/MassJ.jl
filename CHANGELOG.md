# MassJ.jl changelog

## Version `v2.0.0`

A type-system unification. This is a **breaking** release: the public container
types changed shape. Functionality is otherwise unchanged.

* ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The two spectrum
  types `MSscan` (single scan) and `MSscans` (composite) are merged into a single
  type, `MSscans`. Every per-scan provenance field (`num`, `rt`, `level`,
  `precursor`, `polarity`, `activationMethod`, `collisionEnergy`, `chargeState`,
  `driftTime`, `compensationVoltage`) is now a `Vector`. A raw scan is simply an
  `MSscans` with length-1 provenance, and the variance accumulator `s` is empty
  until two or more scans are combined.

* ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The three trace
  types `Chromatogram`, `Mobilogram`, and `Ionogram` are merged into a single
  type, `IonCurrent`, with an `axis` field (`:rt`, `:drift`, or `:cv`) recording
  what the abscissa `x` represents. The maximum ion current is now obtained with
  the function `maxic(trace)` instead of the stored `.maxic` field.

* ![BREAKING][badge-breaking] `load` now returns the same type for every
  multi-spectrum format: an `MSrun` (an `AbstractVector{MSscans}`). Previously
  only mzML returned `MSrun` while mzXML/MGF/MSP/imzML returned a bare
  `Vector{MSscans}`. For the non-mzML formats the `MSrun`'s metadata and
  chromatogram fields are empty. A TXT file (a single spectrum) and a spectrum
  saved as a scalar still load as a bare `MSscans`. The filter/processing
  functions (`extract`, `chromatogram`, `average`, `retention_time`, `smooth`,
  `centroid`, `baseline_correction`) accept any `AbstractVector{MSscans}`, so
  they work on an `MSrun` directly.

* ![Feature][badge-feature] `load_usi`: retrieve a single spectrum by its Universal Spectrum
  Identifier (USI) from a PROXI-compliant server (default the ProteomeCentral aggregator),
  returning an `MSscans` with precursor / charge / activation / collision-energy /
  isolation-window metadata. Lets a spectrum be loaded by reference from a public proteomics
  or metabolomics repository rather than from a local file. Adds `JSON` and `Downloads`
  dependencies.

* ![Feature][badge-feature] New filters for the per-spectrum information exposed by the mzML
  reader, all composable with the existing filters in a single pass: `ChargeState`,
  `SpectrumType` (`:centroid`/`:profile`), `MobilityType`, a generic `MetaData(key[, value])`
  filter over any acquisition cvParam (exact / range / substring / presence), and
  `InstrumentConfig` (match by configuration id, cvParam name, or component type, resolved
  through the run-level instrument table). Each spectrum's `instrumentConfigurationRef` is now
  read into its `metadata`.

* ![Feature][badge-feature] Writers for the remaining open text formats — `save_mgf` (Mascot
  Generic Format), `save_msp` (NIST spectral library), and `save_txt` (two-column m/z/intensity,
  single spectrum) — wired into `save` by extension. Read/write is now symmetric for every
  format except imzML (still read-only).

* ![Feature][badge-feature] Statistical decomposition of chimeric (multiplexed) tandem mass
  spectra, after Truong, Nahon & Giuliani (*J. Am. Soc. Mass Spectrom.* 2026, 37, 649):
  `abundance_matrix` (series → binned abundance matrix + TIC), `partial_correlation` (TIC-
  controlled partial Pearson, via `Statistics`), `cmi_matrix` (conditional mutual information,
  via `EntropyInvariant`), `cluster_ions` (Ward hierarchical clustering + elbow, via
  `Clustering`), and `cluster_spectra` (per-cluster reconstructed spectra). `cmi_matrix` and
  `cluster_ions` ship as **package extensions** — they activate when `EntropyInvariant` /
  `Clustering` are loaded, so neither becomes a hard dependency.

* ![Feature][badge-feature] Continuous Wavelet Transform centroiding (`centroid(scan;
  method = CWT())`), a Ricker-wavelet ridge-line peak detector after Du, Kibbe & Lin
  (*Bioinformatics* 2006). Replaces the former unimplemented `CWT` stub; joins the existing
  `SNRA` and `TBPD` centroiding methods.

* ![Feature][badge-feature] imzML writer (`save_imzml`) for imaging mass spectrometry — emits the
  paired `.imzML` (XML) + `.ibd` (binary, 16-byte UUID + raw little-endian arrays + SHA-1
  checksum) in processed mode, preserving per-pixel coordinates. Read **and** write are now
  supported for every format MassJ reads.

### Migration from v1.x

| v1.x                                   | v2.0                                       |
|----------------------------------------|--------------------------------------------|
| `scan::MSscan`                         | `scan::MSscans` (provenance is now a vector)|
| `scans[1].num` → `1`                   | `scans[1].num` → `[1]` (use `scans[1].num[1]`) |
| `chrom::Chromatogram`, `chrom.rt`      | `chrom::IonCurrent`, `chrom.x`             |
| `chrom.maxic`                          | `maxic(chrom)`                             |
| `Mobilogram`, `Ionogram`               | `IonCurrent` with `axis = :drift` / `:cv`  |

The scalar `MSscans(num, rt, tic, mz, int, level, …)` constructors mirror the
former `MSscan` signatures, so code that *constructed* scans needs only the type
name changed; code that *read* a scalar provenance field must index `[1]`.


## Version `v0.3.1`

* [Enhacement][badge-enhancement]: documentation has been updated.

* ![Feature][badge-feature] `load`: mzxml import now support zlib compression.

* ![Fix][badge-bugfix] Bug fix in mzXML loading: some mzXML contains a precursorMz entry without any ActivationMethod field. A test to ActivationMethod has been added in mzxml.jl.

* ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] Renaming to `MassJ.jl`. This short name should be easier to remember (as it refers to MSn).


## Version `v0.3.0`

* ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] Renaming: `MSJ.jl` to be more inline with the package naming guidelines.

* ![Feature][badge-feature] `formula`: public function. Parse a string containing the chemical formula to a dict{String, Int}.

* ![Feature][badge-feature] `isotopic_distribution`: public function. Takes either the chemical formula or a dict{String,Int} and returns a set of isotopologues.

* ![Feature][badge-feature] `masses`: public function. Takes either the chemical formula or a dict{String,Int} and returns the average, nominal and isotopic masses.
* ![Feature][badge-feature] `simulate`: public function. From an isotpic distribution simulate a mass spectrum based on a resolution and a peak shape.

## Version `v0.2.0`
* ![BREAKING][badge-breaking] `msfilter` function is renamed `average`.

* ![Feature][badge-feature] `baseline_correction`: public function. Corrects the baseline using TopHat, Iterative polynomial smoothing algorithm (IPSA) or Locally weighted error sum of squares regression (LOESS).

* ![Feature][badge-feature] `extract`: public function. Extract a subset of a ms data.


## Version `v0.1.0`
* ![Feature][badge-feature] Supported files: `mzXML`.

* ![Feature][badge-feature] `info`: public function that reads the content of a file, but without loading the data.

* ![Feature][badge-feature] `load`: public function that loads the content of a file.

* ![Feature][badge-feature] `chromatogram`: public function that retrieves chromatograms from the content of a file.

* ![Feature][badge-feature] `msfilter`: public function, returns the average mass spectrum either from a file or from a variable.

* ![Feature][badge-feature] Filtering: msfilter and chromatogram by scan number, MS level, Polarity, Activation Method, Activation Energy, Precursor, retention time, ion current..

* ![Feature][badge-feature] `centroid`: public function, performs peak picking from a file or a variable.

* ![Feature][badge-feature] `smooth`: public function, smooth MS data on variable.




[github-1148]: https://github.com/JuliaDocs/Documenter.jl/issues/1148
[github-1189]: https://github.com/JuliaDocs/Documenter.jl/pull/1189


[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/bugfix-purple.svg
[badge-security]: https://img.shields.io/badge/security-black.svg
[badge-experimental]: https://img.shields.io/badge/experimental-lightgrey.svg
[badge-maintenance]: https://img.shields.io/badge/maintenance-gray.svg

<!--

# Badges

![BREAKING][badge-breaking]
![Deprecation][badge-deprecation]
![Feature][badge-feature]
![Enhancement][badge-enhancement]
![Bugfix][badge-bugfix]
![Security][badge-security]
![Experimental][badge-experimental]
![Maintenance][badge-maintenance]
-->
