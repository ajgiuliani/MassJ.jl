# Energy-resolved yields

The [`yields`](@ref) function builds a [`MassJ.YieldCurve`](@ref) from a *series of
mass-spectrum files* indexed by an external parameter — typically a photon energy
(UVPD, action spectroscopy), a wavelength, or a collision energy (CID breakdown
curves). For each file, the spectrum is averaged once and a set of
[`MassJ.Peak`](@ref) windows are integrated; the result is a matrix of peak
intensities versus that parameter.


## Defining the peak windows

Peaks come in two flavours.

### [`MassJ.Peak`](@ref) — fixed window

A `Peak` carries an explicit m/z window `(mz1, mz2)` plus a label. The same
window is integrated across every spectrum in the series. Three constructors:

```julia
peaks = [
    MassJ.Peak(110.0, 111.0, "fragment_a"),     # explicit bounds (bounds swap if reversed)
    MassJ.Peak(150.5, "fragment_b"; tol = 0.5), # mz ± Δm/z  (absolute)
    MassJ.Peak(500.0, "precursor";  ppm = 5.0), # mz ± mz·ppm·1e-6 (high-res)
]
```

Use this when calibration is stable across the series.

### [`MassJ.TargetPeak`](@ref) — located per file

A `TargetPeak` carries a target m/z plus a search half-width (`tol` or `ppm`).
For each spectrum, the actual peak position is located in `[mz - tol, mz + tol]`
and a window is derived from it according to `method`:

| method        | description                                                          |
|---------------|----------------------------------------------------------------------|
| `:local_max`  | (default) `argmax(int)` in the search window; window is ±`tol` around it |
| `:edges`      | start at the local max, walk outward while `int > edges * peak_max` (default `edges = 0.1`) |
| `:centroid`   | run [`MassJ.centroid`](@ref) on the averaged spectrum and pick the strongest centroid in the search window |

```julia
peaks = [
    MassJ.TargetPeak(110.5, "fragment_a"; tol = 0.5),                       # B
    MassJ.TargetPeak(500.05, "precursor"; ppm = 5.0, method = :edges),      # C
    MassJ.TargetPeak(195.09, "M+H";       tol = 0.5, method = :centroid),   # D
]
```

Use `TargetPeak` when calibration drifts across the series, when peak widths
vary (`:edges`), or when you want the package's centroiding to lock onto each
peak (`:centroid`). The located m/z for each file/peak is stored in
`yc.found_mz` for verification.

For `:centroid`, pass the centroiding parameters to [`yields`](@ref) via
`centroid_method`. The default is `MassJ.SNRA(1.0, 100)`; for high-resolution
data with many close peaks, prefer `MassJ.TBPD(:gauss, R, threshold)` with `R`
matching the instrument resolution.

### Isotope clusters and formula-driven peaks

A `TargetPeak` can bundle **several m/z under one label** — typically the
isotopologues of one species — so [`yields`](@ref) reports a single column for the
whole pattern. The sub-window integrals are summed *inside each scan* before the
error is taken, so the pattern-level uncertainty is correct (the isotopes
co-vary), and overlapping sub-windows are merged so a shared region is counted
once. Three forms:

```julia
# explicit cluster of m/z
MassJ.TargetPeak([442.01, 443.01, 444.01], "Nd cluster"; tol = 0.02)

# from a chemical formula — the isotopologue m/z come from `isotopic_distribution`
# up to cumulative probability `p_target`; with `adduct = ""` the m/z are the
# neutral isotopologue masses divided by |charge|
MassJ.TargetPeak("Nd(NO3)4", "Precursor"; charge = -1, adduct = "",
                 p_target = 0.999, tol = 0.2)

# with an adduct (charge comes from the adduct, m/z via `adduct_mz`)
MassJ.TargetPeak("C6H12O6", "glucose"; adduct = "[M+H]+", tol = 0.01)
```

For a cluster, `method` is `:fixed` (default — fixed windows at the theoretical
m/z) or `:anchor` (locate the **anchor** isotopologue in each spectrum and shift
the whole pattern by the calibration offset, preserving the isotope spacing;
falls back to fixed when the anchor is absent). The anchor defaults to the
most-abundant isotopologue for a formula (`anchor = :max` — robust for heavy ions
whose monoisotopic peak is weak); use `anchor = :mono` for the lowest m/z, or pass
an explicit m/z. Because the whole pattern is shifted rigidly by one offset,
isotopologues lighter than the anchor reconstitute correctly below it.

### Loading peaks from CSV

[`MassJ.read_peaklist`](@ref) auto-detects three CSV layouts by column count:

| Cols | Header example                  | Result                                                         |
|------|---------------------------------|----------------------------------------------------------------|
| 2    | `mz,label`                      | [`TargetPeak`](@ref) using the `tol`/`ppm`/`method` kwargs     |
| 3    | `mz1,mz2,label`                 | [`Peak`](@ref) (fixed window — legacy form)                    |
| 4    | `mz,tol,method,label`           | [`TargetPeak`](@ref) with per-row tolerance and method         |

```julia
# 3-column legacy form
peaks = MassJ.read_peaklist("peaks.csv")

# 2-column form: targets only, with global defaults
peaks = MassJ.read_peaklist("targets.csv"; tol = 0.5, method = :local_max)

# 4-column form: per-row method
peaks = MassJ.read_peaklist("targets_full.csv")
```

The header row is optional and detected automatically. Example 4-column file:
```
mz,tol,method,label
110.5,0.5,local_max,fragment_a
500.05,0.005,edges,precursor
```


## Building a YieldCurve

The [`yields`](@ref) function has two methods.

**From an explicit list of files** — when you control the x-axis exactly:
```julia
files = ["scan_3.0eV.mzML", "scan_3.5eV.mzML", "scan_4.0eV.mzML"]
yc = yields(files, peaks;
            x      = [3.0, 3.5, 4.0],
            xlabel = "photon energy (eV)")
```

**From a directory** — when files are named on a regular grid and `x = x0 + step·i`:
```julia
yc = yields("data/UVPD/", peaks;
            x0     = 3.0,
            step   = 0.1,
            xlabel = "photon energy (eV)")
```
Files in the directory are picked in *natural-sort* order (so `scan2.mzML` comes
before `scan10.mzML`). Supported extensions are mzXML, mzML, MGF, MSP, imzML, and
TXT.

### Filtering and preprocessing each file

Any [`MassJ.FilterType`](@ref) passed after `peaks` is applied to each file's scans
before averaging and integration, so the yield is built from only the scans you
want (one MS level, a retention-time window, an activation method, …):

```julia
yc = yields(files, peaks, Level(1), RT(5.0, 15.0); x = energies)
yc = yields("data/UVPD/", peaks, Level(1); x0 = 3.0, step = 0.1)
```

For full control over per-point preprocessing (filtering, smoothing, baseline
correction, calibration, …) build the scan lists yourself and pass a *series* —
one `Vector{MSscans}` per point:

```julia
series = [baseline_correction.(smooth.(load(f))) for f in files]
yc = yields(series, peaks; x = energies, labels = string.(energies))
```

Both forms preserve the per-scan error model described in the
**Uncertainties** section below.

The result is a [`MassJ.YieldCurve`](@ref) holding the energy axis `yc.x`, the
matrix `yc.yields[file, peak]`, the per-file `yc.tic` (sum across peak windows),
the matrix `yc.found_mz[file, peak]` (`NaN` for fixed `Peak`s, the located m/z
for each `TargetPeak`), and the labels and nominal windows used. Two parallel
matrices, `yc.yields_err` and `yc.tic_err`, carry the propagated 1-σ
uncertainties — see [Uncertainties](@ref) below.

When any peak in the list is a `TargetPeak` with `method = :centroid`,
[`yields`](@ref) also accepts a `centroid_method` keyword forwarded to
[`MassJ.centroid`](@ref):

```julia
yc = yields("data/", peaks; x0 = 3.0, step = 0.1,
            centroid_method = MassJ.TBPD(:gauss, 4500.0, 0.2))
```


## Normalization

Normalization steps are *post-processing* functions that return a new
`YieldCurve` (the input is not mutated). They can be composed:

```julia
yc_norm = yc |> normalize_tic |> y -> normalize_external(y, "flux.txt")
```

[`normalize_tic`](@ref) divides each row of `yc.yields` by that row's TIC so the
peaks sum to 1 per energy step (relative branching ratios). The raw TIC column
is preserved.

[`normalize_external`](@ref) divides each row (peaks and TIC) by an **external
reference** — a quantity that is *not* contained in the spectra — interpolated to
that row's `x` value from a text file. The reference can be anything the yields
should be normalised against: photon flux, laser power or pulse energy, ion-source
current, detector efficiency, transmission, or acquisition time. See the
[External-reference file format](@ref) section below.

```julia
yc_n = normalize_external(yc_tic, "flux.txt")                    # 10% default σ
yc_n = normalize_external(yc_tic, "power.txt"; err_pct = 0.05)   # override default
yc_n = normalize_external(yc_tic, "ref.txt";   skipstart = 3)    # force header skip
yc_n = normalize_external(yc_tic, "ref.txt";   extrapolate = :line)  # linear extrap.
```

By default (`extrapolate = :clamp`), `x` values outside the reference file's range
are clamped to the nearest endpoint and a warning is emitted. Pass
`extrapolate = :line` to linearly extrapolate using the slope of the nearest
segment (the same behaviour as `Interpolations.LinearInterpolation(...; extrapolation_bc = Line())`).
The same scheme is applied to the σ column, with `abs(...)` to keep uncertainties
non-negative.

!!! note "Renamed from `normalize_flux`"
    This function was previously `normalize_flux` (with a `flux_err_pct` keyword).
    The old name is kept as a deprecated alias that forwards to
    `normalize_external`, so existing scripts keep working.

Each normalization records itself in `yc.metadata` so the provenance of a curve
is preserved.


## External-reference file format

[`normalize_external`](@ref) is liberal in what it accepts. Three layouts are
recognised, distinguished by the column count (`v` is the reference value):

| Cols | Meaning                  | Where σ comes from                                |
|------|--------------------------|---------------------------------------------------|
| 2    | `x  v`                   | `err_pct * v` (default 10%)                       |
| 3    | `x  v  σ_v`              | the third column (per-row)                        |
| ≥3   | `x  v  <non-numeric>`    | falls back to `err_pct * v` for that row          |

Parsing rules:

* Lines starting with `#` are treated as comments and stripped anywhere in the
  file (header block, between data rows, or as a trailing `# ...` on a data
  line).
* Leading non-numeric rows (text column headers like `energy flux`) are
  auto-detected and skipped.
* Whitespace separates columns (any mix of spaces and tabs).
* When a 3rd column is present but contains non-numeric text on a given row —
  for example a timestamp on each line of a DESIRS beamline flux log — that
  row's σ silently falls back to `err_pct * v`. No need to pre-process
  the file.
* Pass `skipstart = N` to discard `N` *physical* lines from the top before
  parsing, for the rare cases where a numeric-looking row at the top is not
  data (e.g. a unit/scale row `1 1`).

Example DESIRS-style file that just works:
```
####################################################
# File created by DESIRS on Fri Apr 17 12:27:46 2026
# Scan parameters …
####################################################
Energy (eV)    K6514 (A)    Time
3.999993    2.169652e-010    Fri Apr 17 12:27:54 2026
4.099999    3.912680e-010    Fri Apr 17 12:28:02 2026
4.200016    7.635849e-010    Fri Apr 17 12:28:10 2026
```
The `####` and `#`-prefixed comment block is stripped, the column header line
is auto-detected and skipped, and the timestamp column on each data line is
recognised as non-numeric → σ_φ defaults to 10% of φ.


## Uncertainties

For each peak window, [`yields`](@ref) integrates the window **in every
individual scan** and takes the standard error of the mean across the
per-scan areas:

```
area_k    = ∫ scan_k.int dm/z   over the window     (k = 1 .. N)
yield     = mean(area_k)
yields_err = std(area_k) / sqrt(N)
```

The per-file peak window is fixed (located once on the averaged spectrum so
the position is stable across scans), but the integration runs on every
underlying scan from `load(file)`. The resulting `yields_err` captures the
true scan-to-scan variability of the peak area — peak-height fluctuations,
shape changes, baseline drift, etc. — without the optimistic
bin-independence assumption that diagonal-only error propagation makes.

The combined uncertainty on each row's total is `yc.tic_err[i] = sqrt(Σ_p
yields_err[i, p]²)`. For a single-scan input (`MSscans` or any file
containing only one scan), there is no scan-to-scan variation to estimate
and the corresponding entries are `NaN`.

!!! note "Why per-scan, not per-bin?"
    The natural-looking alternative — derive a per-m/z-bin standard error
    from the Welford accumulator stored in `MSscans.s` and propagate it
    through the trapezoidal integral — assumes the m/z bins fluctuate
    *independently* across scans. They don't: across a real peak, adjacent
    bins rise and fall together. The diagonal-only formula then drops the
    (large positive) off-diagonal covariance terms and **underestimates the
    true uncertainty by an order of magnitude or more**. Per-scan area
    integration captures the full per-peak covariance without making any
    independence assumption.

The errors propagate through subsequent normalization:

* [`normalize_tic`](@ref) — `σ(y/T) = (y/T)·sqrt((σ_y/y)² + (σ_T/T)²)` (the
  correlation between `y` and `T = Σy` is ignored, the usual first-order
  approximation).
* [`normalize_external`](@ref) — `σ(y/v) = (y/v)·sqrt((σ_y/y)² + (σ_v/v)²)`. The
  reference uncertainty `σ_v` is taken from a third column in the reference file
  when present; otherwise it is `err_pct * v` (default 10%). Use the keyword
  `err_pct = 0.05` to override:
  ```julia
  yc_n = normalize_external(yc_tic, "ref.txt"; err_pct = 0.05)
  ```

Plotting picks up the ribbon automatically — see [Plotting](#Plotting).


## Selecting peaks for a plot

[`drop_peaks`](@ref) returns a new `YieldCurve` with one or more peaks removed,
keeping everything else (including the `tic` reference and `found_mz`)
identical. This is the usual workflow when one peak — typically the precursor —
swamps the axes:

```julia
plot(drop_peaks(yc, "precursor"))
plot(drop_peaks(yc, ["precursor", "solvent"]))
```

Labels that don't match any peak in `yc.labels` are silently ignored. The
`tic` field is left unchanged so the dropped curve still references the
original total — useful if you also want to normalize.


## Plotting

A [Plots.jl](https://github.com/JuliaPlots/Plots.jl) recipe is provided — one
line per peak, x-axis from `yc.x`, legend from `yc.labels`. When
`yc.yields_err` contains finite values, 1-σ ribbons are drawn around each line
with `fillalpha = 0.15`. Override either by passing the same keywords to `plot`.

```julia
using Plots
plot(yc_norm)                                 # ribbon by default
plot(yc_norm; ribbon = nothing)               # no ribbon
plot(yc_norm; fillalpha = 0.30)               # darker ribbon
```


## Interoperability: DataFrames and Tables.jl

A [`MassJ.YieldCurve`](@ref) is a [Tables.jl](https://github.com/JuliaData/Tables.jl)
source. The interface is provided by a package extension that loads automatically
as soon as `Tables` — or anything built on it, such as `DataFrames` or `CSV` — is
in scope. The table has one row per x value, with columns `x`, one per peak
label, a matching `<label>_err`, then `TIC`, `TIC_err`, and `file`:

```julia
using DataFrames
df = DataFrame(yc)                  # straight into a DataFrame

using CSV
CSV.write("yields.csv", yc)         # any Tables sink
```

This lets a yield curve flow directly into the rest of the Julia data ecosystem
(DataFrames, query/split-apply-combine, Arrow, StatsPlots …).


## Writing to CSV

[`MassJ.write_csv`](@ref) writes a compact CSV with header
`<xlabel>, <peak labels…>, TIC` (no error columns) — independent of the Tables
interface above; use whichever layout suits your downstream workflow:

```julia
MassJ.write_csv(yc_norm, "yields.csv")
```


## Command-line interface

A thin CLI wrapper is provided at `scripts/yields.jl` for batch processing
without writing Julia code:

```bash
julia --project=@. scripts/yields.jl \
      --input data/UVPD/ \
      --list  peaks.csv \
      --energy 3.0 --step 0.1 \
      --xlabel "photon energy (eV)" \
      --tic \
      --flux flux.txt \
      --output yields.csv
```

Run `julia scripts/yields.jl --help` for the full list of flags.


## Low-level primitive

[`integrate_window`](@ref) is the building block used internally by `yields` and
can also be called directly on any [`MassJ.MScontainer`](@ref):

```julia
spec = average("scan.mzML")
area = integrate_window(spec, 150.0, 152.0)
```

It performs trapezoidal integration over the window and returns `0.0` when fewer
than two sample points fall inside.
