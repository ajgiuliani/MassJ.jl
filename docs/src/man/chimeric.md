# Chimeric tandem mass spectra

When several precursor ions are activated simultaneously (for example by photon
activation over a wide isolation window), the resulting tandem mass spectrum is
*chimeric* — it mixes fragments from every precursor. MassJ implements the
data-driven framework of Truong, Nahon & Giuliani
([*J. Am. Soc. Mass Spectrom.* **2026**, 37, 649](https://doi.org/10.1021/jasms.5c00358))
to statistically resolve such spectra into their per-precursor contributions,
without prior knowledge of the precursor identities.

The idea: record many tandem mass spectra while the activation conditions
fluctuate, so that fragments from the *same* precursor vary together across
acquisitions while fragments from *different* precursors do not. Scoring the
co-variation of every pair of ions, then clustering, recovers the groups.

## Workflow

```julia
using MassJ
using EntropyInvariant      # enables mi / cmi     (package extension)
using Clustering            # enables cluster_ions (package extension)

scans = load("chimeric_series.mzML")          # N acquisitions (an MSrun)

# 1. Reduce the series to an abundance matrix (threshold + bin + align)
am = abundance_matrix(scans; binsize = 1.0, threshold = 0.01)
#   am.matrix :: N×M   am.mz :: M bin centres   am.tic :: N per-acquisition TIC

# 2. Score every ion pair (pure functions on the matrix — binning already done)
P = pcor(am.matrix, am.tic)          # partial correlation, conditioned on the TIC
#   cor(am.matrix)                   # plain correlation
#   cmi(am.matrix, am.tic)           # conditional mutual information

# 3. Cluster the ions and reconstruct one spectrum per precursor
labels  = cluster_ions(P; kind = :correlation) # Ward linkage, elbow chooses k
spectra = cluster_spectra(am.mz, am.matrix, labels)   # Vector{MSscans}

using MassJ.plots, Plots
plot(spectra[1])                                # the first resolved precursor
```

## The abundance matrix

[`abundance_matrix`](@ref) turns the series into the `N × M` matrix that the
statistics operate on. Each spectrum is thresholded (peaks below `threshold` ×
its base peak are dropped, default 1 %) and its surviving intensities are summed
into `binsize`-wide m/z bins on a common axis. Columns are the *variables*
``X_i`` (one binned m/z), rows are the acquisitions. The per-acquisition total
``I = \sum_i X_i`` is returned as `am.tic` and is used as the control /
conditioning variable. Empty bins (no varying signal) are dropped by default.

## Choosing the metric

Five score functions cover the useful metrics. They are pure — matrix in, `M × M` score
matrix out — and the conditioned ones take the conditioning variable explicitly:

| | unconditioned | conditioned |
|---|---|---|
| covariance | `cov(X)` | *`pcov` — patent-encumbered* |
| correlation | `cor(X)` | [`pcor`](@ref)`(X, control)` |
| mutual information | [`mi`](@ref)`(X)` | [`cmi`](@ref)`(X, control)` |

`cov`/`cor` are reused from `Statistics`; [`pcor`](@ref) is covariance-based (linear,
signed); [`mi`](@ref)/[`cmi`](@ref) come from `EntropyInvariant` (nonlinear, always ≥ 0,
sign-blind). Conditioning on `control = am.tic` removes the intensity fluctuation common
to all ions, leaving a positive score for ions of the same precursor and a negative one
between precursors (for the signed measures). Correlation is ideal for clean
correlated/anticorrelated cases (two precursors); mutual information detects dependence
of any shape, useful when a simple split breaks down (three or more precursors). The
`control` can be any per-acquisition vector — an internal-standard intensity, an
instrument parameter — not just the TIC. Conditioned covariance (`pcov`) is absent
because it is patent-encumbered.

## Clustering and reconstruction

[`cluster_ions`](@ref) applies Ward hierarchical agglomerative clustering to the
score matrix (converting it to a dissimilarity according to `kind`) and returns a
cluster label per ion. The number of clusters is chosen automatically by the
elbow of the within-cluster sum of squares, or fixed with `nclusters`.
[`cluster_spectra`](@ref) then averages each cluster's bins across the
acquisitions to produce one [`MassJ.MSscans`](@ref) per precursor.

## Significance and hands-off correspondence

[`jackknife`](@ref) and [`permutation`](@ref) attach a *per-edge p-value* to a score, so
the grouping can be gated by statistical confidence rather than a hand-picked threshold.
Both are higher-order — you pass the score function — and return `(estimate, pvalue)`:

- [`jackknife`](@ref)`(score, X[, control])` — a delete-one jackknife: a small p marks a
  co-variation supported by the whole series rather than a single outlier acquisition
  (the failure mode that inflates raw covariance). Valid for the covariance/correlation
  family.
- [`permutation`](@ref)`(score, X[, control]; nperm)` — a column-permutation p-value, the
  route for `mi`/`cmi`, whose null is not zero-centred.

[`fdr_adjust`](@ref) turns the p-value matrix into a Benjamini-Hochberg gate at a chosen
FDR level, keeping only edges that survive multiple-comparison correction over the
``M(M-1)/2`` pairs. [`correspondence`](@ref) chains all of it into one call — abundance
matrix → score → significance → FDR gate → clustering → per-precursor spectra:

```julia
res = correspondence(scans; measure = :pcor, control = :tic, method = :jackknife, q = 0.05)
res.labels        # cluster label per ion (as from cluster_ions)
res.spectra       # Vector{MSscans}, one resolved precursor each
res.significant   # BitMatrix of the FDR-surviving edges
```

`measure` picks the score (`:cov`, `:cor`, `:pcor`, `:mi`, `:cmi`); `control` is the
conditioning variable for `:pcor`/`:cmi` (`:tic`, `:none`, or a vector); `:mi`/`:cmi`
need `method = :permutation`. If no edge survives the gate it warns and falls back to the
ungated score, so you always get a result.

## Covariance maps and fragment pairs

Instead of *grouping* ions, the covariance-mapping route reports individual correlated
*pairs*. `cov` builds the map and [`covariance_features`](@ref)
picks its peaks and quantifies each:

```julia
feats = covariance_features(scans; normalize = :tic, nfeatures = 50)
# each MapFeature :: (mz1, mz2, value, volume, snr), sorted by jackknife S/N
```

`normalize = :tic` divides each spectrum by its total ion current before the
covariance — a simple, unconditioned way to damp the common intensity fluctuation (it
uses no conditioning model). Internally [`map_features`](@ref) does the 2-D
peak-picking after [`cut_autocorrelation`](@ref) suppresses the diagonal; each pair's
volume is a windowed 2-D integral of the covariance, and its `snr` is a leave-one-out
jackknife of that volume. `map_features` works on any symmetric map, so you can also
pick pairs off a `pcor` or `cmi` map. Visualise the map with the
classic marginal-spectrum figure [`MassJ.plots.association_map_marginal`](@ref):

```julia
using MassJ.plots, Plots
am = abundance_matrix(scans; binsize = 1.0)
A  = cov(am.matrix ./ max.(am.tic, 1))    # TIC-normalised
association_map_marginal(A, am.mz, vec(sum(am.matrix, dims = 1)))
```

## Optional dependencies

The scoring steps [`mi`](@ref)/[`cmi`](@ref) and the clustering step
[`cluster_ions`](@ref) are provided as **package extensions**: they become
available only once `EntropyInvariant` and `Clustering` are loaded in your
session. `abundance_matrix`, `cov`/`cor`/`pcor`, and `cluster_spectra` are
always available (they depend only on `Statistics`). Calling an extension
function before loading its backing package raises an error telling you which
`using` statement to add.
