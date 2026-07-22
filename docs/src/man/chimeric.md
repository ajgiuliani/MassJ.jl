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
using EntropyInvariant      # enables cmi_matrix   (package extension)
using Clustering            # enables cluster_ions (package extension)

scans = load("chimeric_series.mzML")          # N acquisitions (an MSrun)

# 1. Reduce the series to an abundance matrix (threshold + bin + align)
am = abundance_matrix(scans; binsize = 1.0, threshold = 0.01)
#   am.matrix :: N×M   am.mz :: M bin centres   am.tic :: N per-acquisition TIC

# 2a. Score by partial Pearson correlation (TIC as the control variable)
P = partial_correlation(am.matrix)             # M×M, signed

# 2b. …or by conditional mutual information (nonlinear, sign-blind, ≥ 0)
C = cmi_matrix(am.matrix)                       # M×M, conditioned on the TIC

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

## Choosing the score

| Score | Function | Range | Captures | Best for |
|-------|----------|-------|----------|----------|
| Partial Pearson correlation | [`partial_correlation`](@ref) | −1 … 1 | linear, **signed** co-variation | clean correlated/anticorrelated cases (e.g. two precursors) |
| Conditional mutual information | [`cmi_matrix`](@ref) | ≥ 0 | any (nonlinear) dependence | complex mixtures where sign-based separation breaks down (e.g. three precursors) |

Partial correlation conditions on the TIC to remove the intensity fluctuation
common to all ions, leaving positive scores for ions of the same precursor and
negative scores between precursors. Conditional mutual information likewise
conditions on the TIC but detects relationships of any shape; because it is
sign-blind it is most useful when a simple correlated/anticorrelated split is no
longer adequate.

## Clustering and reconstruction

[`cluster_ions`](@ref) applies Ward hierarchical agglomerative clustering to the
score matrix (converting it to a dissimilarity according to `kind`) and returns a
cluster label per ion. The number of clusters is chosen automatically by the
elbow of the within-cluster sum of squares, or fixed with `nclusters`.
[`cluster_spectra`](@ref) then averages each cluster's bins across the
acquisitions to produce one [`MassJ.MSscans`](@ref) per precursor.

## Significance and hands-off correspondence

Any association can carry a *per-edge significance*, so the grouping can be gated by
statistical confidence rather than a hand-picked score threshold. Two resampling
routes are provided:

- [`jackknife_significance`](@ref) — a delete-one jackknife giving each pair an
  `(estimate, se, z)`. A large `|z|` marks a co-variation supported by the whole
  series rather than a single outlier acquisition (the failure mode that inflates raw
  covariance). Use it for the covariance / correlation family (`:covariance`,
  `:correlation`, `:pcorr`).
- [`permutation_significance`](@ref) — a column-permutation p-value, the correct route
  for a statistic whose null is not zero-centred, notably conditional mutual
  information (`:cmi`).

[`fdr_adjust`](@ref) turns a p-value matrix into a Benjamini-Hochberg gate at a chosen
FDR level, keeping only edges that survive multiple-comparison correction over the
``M(M-1)/2`` pairs. [`correspondence`](@ref) chains all of it into one call —
abundance matrix → score → significance → FDR gate → clustering → per-precursor
spectra:

```julia
res = correspondence(scans; score = :pcorr, significance = :jackknife, q = 0.05)
res.labels        # cluster label per ion (as from cluster_ions)
res.spectra       # Vector{MSscans}, one resolved precursor each
res.significant   # BitMatrix of the FDR-surviving edges
```

With `score = :cmi` it uses the permutation route automatically. If no edge survives
the gate it warns and falls back to the ungated score, so you always get a result.

## Covariance maps and fragment pairs

Instead of *grouping* ions, the covariance-mapping route reports individual correlated
*pairs*. [`covariance_matrix`](@ref) builds the map and [`covariance_features`](@ref)
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
pick pairs off a `partial_correlation` or `cmi_matrix`. Visualise the map with the
classic marginal-spectrum figure [`MassJ.plots.covmap_marginal`](@ref):

```julia
using MassJ.plots, Plots
am = abundance_matrix(scans; binsize = 1.0)
A  = covariance_matrix(am.matrix ./ max.(am.tic, 1))    # TIC-normalised
covmap_marginal(A, am.mz, vec(sum(am.matrix, dims = 1)))
```

## Optional dependencies

The scoring step [`cmi_matrix`](@ref) and the clustering step
[`cluster_ions`](@ref) are provided as **package extensions**: they become
available only once `EntropyInvariant` and `Clustering` are loaded in your
session. `abundance_matrix`, `partial_correlation`, and `cluster_spectra` are
always available (they depend only on `Statistics`). Calling an extension
function before loading its backing package raises an error telling you which
`using` statement to add.
