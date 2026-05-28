# Uncertainties and errors

MassJ derives intensity uncertainties from **replicate scans**: a single raw
spectrum carries no error, but as soon as two or more scans are combined the
package tracks the per-m/z spread. This page shows how to retrieve those errors
after **averaging**, after **spectral arithmetic**, and from an energy-resolved
**yield** curve.

## The uncertainty field and its accessors

Every [`MassJ.MSscans`](@ref) has an uncertainty field `s`. You normally never
read it directly — its meaning depends on the stage of the computation (the
incremental Welford *M2* accumulator while scans are being combined, then the
per-m/z sample standard deviation once finalized by [`average`](@ref)). Use the
accessors instead:

| accessor | returns |
|----------|---------|
| [`nscans`](@ref)  | number of scans averaged into the spectrum (`1` for a raw scan) |
| [`stdev`](@ref)   | per-m/z **sample standard deviation** σ — the spread of the replicate intensities |
| [`sem`](@ref)     | per-m/z **standard error of the mean** σ/√N — the 1-σ uncertainty on each averaged intensity |

All three return an empty vector when no replicate statistics are available
(a single scan, or `average(...; stats = false)`).

## Averaging

[`average`](@ref) combines spectra with an incremental Welford algorithm. With
`stats = true` (the default) the returned spectrum carries the per-m/z standard
deviation, from which the accessors derive everything:

```julia
using MassJ

scans = load("run.mzML")            # several scans
avg   = average(scans, Level(1))    # stats = true by default

n  = nscans(avg)                    # how many scans were averaged
σ  = stdev(avg)                     # per-m/z standard deviation  (length == length(avg.mz))
se = sem(avg)                       # per-m/z standard error of the mean  (σ ./ √n)

# the averaged intensity at point i is  avg.int[i] ± se[i]   (1-σ)
```

Pass `stats = false` to skip the variance computation when you do not need it
(slightly faster); `stdev(avg)` and `sem(avg)` are then empty.

### Plotting the band

The mass-spectrum plot recipe draws a ±1-σ ribbon on request:

```julia
using Plots
plot(avg; band = :sem)                      # ribbon = standard error of the mean
plot(avg; band = :std)                      # ribbon = sample standard deviation
plot(avg; band = :sem, method = :absolute)  # works in absolute mode too
```

`band = :none` (the default) draws no ribbon, and the band is silently ignored
for a single scan that carries no replicate statistics.

## Summing and spectral arithmetic

The overloaded operators `+`, `-`, `*`, `/` (see [Base overloaded operators](@ref))
carry the uncertainty vector through the operation, so the same accessors work on
the result:

```julia
a = average(group_a)        # σ_a stored in a
b = average(group_b)        # σ_b stored in b

s    = a + b                # uncertainty vectors combined
half = a / 2                # uncertainty scaled with the intensity
stdev(s); sem(half)         # retrieved the usual way
```

Mechanically:

* `*` and `/` by a scalar scale the stored uncertainty linearly
  (`a * 2`, `a / 100`).
* `+` and `-` add the two stored uncertainty vectors elementwise, interpolating
  onto a common m/z axis first when the two spectra differ in length.

!!! note "Arithmetic needs replicate statistics to carry an error"
    Two *raw* single scans have an empty `s`, so `stdev(scan1 + scan2)` is empty
    — there is no scan-to-scan information to propagate. When you need an
    uncertainty on a combined spectrum, combine spectra that came from
    [`average`](@ref). For rigorous replicate statistics, prefer building the
    average directly (`average(scans)`) over hand-summing scans.

## Yields

For an energy-resolved [`MassJ.YieldCurve`](@ref), the uncertainty is computed
independently — from the spread of the **per-scan** peak areas — and stored in
two parallel matrices:

```julia
yc = yields("data/UVPD/", peaks; x0 = 3.0, step = 0.1)

yc.yields_err        # nfiles × npeaks : 1-σ on each yield
yc.tic_err           # per-file 1-σ on the row total  (= sqrt(Σ_p yields_err²))
plot(yc)             # 1-σ ribbons drawn automatically
```

The full error model — why it integrates per scan rather than per m/z bin, and
how the errors propagate through [`normalize_tic`](@ref) and
[`normalize_flux`](@ref) — is described on the
[Energy-resolved yields](@ref) page.
