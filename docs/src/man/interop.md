# Interoperability

MassJ's typed containers ([`MassJ.MSscans`](@ref), [`MassJ.MSrun`](@ref)) do not
get in the way of the wider Julia ecosystem. Every machine-learning and statistics
library speaks the same lingua franca — a **feature matrix** of `n_samples ×
n_features` — and the only real obstacle is that spectra live on *heterogeneous
m/z axes*. MassJ bridges that gap explicitly: featurise a spectrum series to a
matrix, run any external algorithm, and map the result straight back. The matrix
never replaces your spectra, so a per-row result (an outlier score, a cluster
label) indexes back into the original series by position.

This mirrors the opt-in extension philosophy used for Tables/Measurements/Unitful:
the typed core stays, and thin converters reach out to matrices and tables.

## Featurising a series

[`featurize`](@ref) bins each spectrum onto a common m/z axis and returns an `N×M`
matrix (one **row per spectrum**, row-aligned to the input) plus the bin centres:

```julia
scans  = load("series/")            # Vector{MSscans}, N spectra
X, mz  = featurize(scans; binsize = 1.0)
size(X)        # (N, M)  — rows = spectra, columns = m/z bins
```

The defaults keep all intensity and all bins, so `X` is a faithful binned image of
the spectra that [`from_matrix`](@ref) can rebuild. `X` drops straight into
DecisionTree.jl, MLJ, IsolationForest, PCA, MassJ's own Clustering.jl extension
([`cluster_spectra`](@ref)), and so on.

## Mapping a result back

Because `X` is row-aligned to `scans`, any per-row result is mapped back with
[`select_spectra`](@ref) (a boolean mask or integer indices):

```julia
scores  = my_outlier_model(X)            # one score per row
inliers = select_spectra(scans, scores .< threshold)   # → Vector{MSscans}
clean   = average(inliers)
```

(`select_spectra`, not `select`, to avoid clashing with `DataFrames.select`.)

## Rebuilding spectra from a matrix

When you transform the intensity matrix in a way [`select_spectra`](@ref) cannot
express — scaling, denoising, imputation, an autoencoder reconstruction —
[`from_matrix`](@ref) turns it back into spectra. Each **row** becomes one
[`MassJ.MSscans`](@ref) on the **binned** m/z axis `mz`; pass `template` to carry
each spectrum's metadata (retention time, scan number, polarity, …) onto the
rebuilt one:

```julia
X, mz   = featurize(scans)
X2      = denoise(X)
rebuilt = from_matrix(X2, mz; template = scans)   # Vector{MSscans} on the binned axis
```

The result lives on the binned axis, so it is a *resampled* view of the originals,
not a bit-identical copy.

## A Tables.jl view

[`spectra_table`](@ref) presents a spectrum series as a
[Tables.jl](https://github.com/JuliaData/Tables.jl)-ready `NamedTuple` — one row
per spectrum, with metadata columns `num`, `rt`, `tic` followed by one column per
m/z bin. It needs no extra dependency on MassJ's side (it mirrors
[`isotope_table`](@ref)):

```julia
using DataFrames
df = DataFrame(spectra_table(scans; binsize = 1.0))
```

For a large number of bins, prefer the raw matrix from [`featurize`](@ref) — a
`NamedTuple` with thousands of fields is costly to build.

Full signatures are on the [References](@ref) page.
