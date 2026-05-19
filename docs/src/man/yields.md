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

The result is a [`MassJ.YieldCurve`](@ref) holding the energy axis `yc.x`, the
matrix `yc.yields[file, peak]`, the per-file `yc.tic` (sum across peak windows),
the matrix `yc.found_mz[file, peak]` (`NaN` for fixed `Peak`s, the located m/z
for each `TargetPeak`), and the labels and nominal windows used.

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
yc_norm = yc |> normalize_tic |> y -> normalize_flux(y, "flux.txt")
```

[`normalize_tic`](@ref) divides each row of `yc.yields` by that row's TIC so the
peaks sum to 1 per energy step (relative branching ratios). The raw TIC column
is preserved.

[`normalize_flux`](@ref) divides each row (peaks and TIC) by the photon flux at
that row's `x` value, linearly interpolated from a two-column text file
(2 header lines, then `energy  flux`). Out-of-range `x` values are clamped to
the nearest flux and a warning is emitted.

Each normalization records itself in `yc.metadata` so the provenance of a curve
is preserved.


## Plotting

A [Plots.jl](https://github.com/JuliaPlots/Plots.jl) recipe is provided — one
line per peak, x-axis from `yc.x`, legend from `yc.labels`:

```julia
using Plots
plot(yc_norm)
```


## Writing to CSV

[`MassJ.write_csv`](@ref) writes a CSV with header
`<xlabel>, <peak labels…>, TIC`:

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
