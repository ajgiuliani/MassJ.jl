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

* ![BREAKING][badge-breaking] `load` on a multi-spectrum mzML file always returns
  an `MSrun` (which is an `AbstractVector{MSscans}`); a single spectrum saved as a
  scalar still round-trips back to a bare `MSscans`.

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
