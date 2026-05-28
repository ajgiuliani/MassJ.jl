"""
Statistical decomposition of chimeric (multiplexed) tandem mass spectra.

Implements the data-driven workflow of Truong, Nahon & Giuliani
(*J. Am. Soc. Mass Spectrom.* **2026**, 37, 649): a series of `N` tandem mass
spectra recorded during photon activation is reduced to an abundance matrix
(rows = acquisitions, columns = binned m/z variables); the columns are scored
against one another by partial Pearson correlation ([`partial_correlation`](@ref))
or conditional mutual information ([`cmi_matrix`](@ref)); hierarchical clustering
of the score matrix ([`cluster_ions`](@ref)) groups the ions by precursor, and
[`cluster_spectra`](@ref) reconstructs one mass spectrum per cluster.

[`cmi_matrix`](@ref) requires the `EntropyInvariant` package and
[`cluster_ions`](@ref) the `Clustering` package; both are loaded on demand as
package extensions (`using EntropyInvariant` / `using Clustering`).
"""

export abundance_matrix, partial_correlation, cmi_matrix, cluster_ions, cluster_spectra


"""
    abundance_matrix(spectra; binsize = 1.0, threshold = 0.01, mzrange = nothing,
                     drop_constant = true)

Reduce a series of spectra to the abundance matrix used for the statistical
analysis. `spectra` is any `AbstractVector{MSscans}` (e.g. an [`MSrun`](@ref)).

Each spectrum is thresholded (intensities below `threshold` × its base peak are
discarded) and its surviving peaks are summed into `binsize`-wide m/z bins on a
common axis spanning `mzrange` (or the global m/z range when `nothing`). The
result is returned as a named tuple:

* `matrix` — `N × M` `Matrix{Float64}`, one row per acquisition, one column per
  m/z bin (the statistical *variables* ``X_i``).
* `mz`     — length-`M` vector of bin centres.
* `tic`    — length-`N` per-acquisition total ion current ``I = \\sum_i X_i``
  (the control / conditioning variable).

When `drop_constant = true` (default) bins that carry no varying signal
(zero variance across acquisitions) are removed, since they cannot be correlated.
"""
function abundance_matrix(spectra::AbstractVector{MSscans};
                          binsize::Real = 1.0,
                          threshold::Real = 0.01,
                          mzrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
                          drop_constant::Bool = true)
    N = length(spectra)
    N == 0 && error("abundance_matrix: no spectra supplied.")
    bs = float(binsize)

    if mzrange === nothing
        nonempty = (s for s in spectra if !isempty(s.mz))
        lo = minimum(minimum(s.mz) for s in nonempty)
        hi = maximum(maximum(s.mz) for s in nonempty)
    else
        lo, hi = float(mzrange[1]), float(mzrange[2])
    end

    M = max(1, Int(cld(hi - lo, bs)))
    centers = [lo + (b - 0.5) * bs for b in 1:M]
    X = zeros(Float64, N, M)

    for (n, s) in enumerate(spectra)
        isempty(s.int) && continue
        thr = threshold * maximum(s.int)
        @inbounds for k in eachindex(s.mz)
            v = s.int[k]
            v >= thr || continue
            (lo <= s.mz[k] <= hi) || continue
            b = clamp(Int(fld(s.mz[k] - lo, bs)) + 1, 1, M)
            X[n, b] += v
        end
    end

    if drop_constant
        keep = [var(@view X[:, j]) > 0 for j in 1:M]
        X = X[:, keep]
        centers = centers[keep]
    end

    tic = vec(sum(X, dims = 2))
    return (matrix = X, mz = centers, tic = tic)
end


"""
    partial_correlation(X; control = vec(sum(X, dims = 2)))

Symmetric `M × M` matrix of partial Pearson correlations between the columns of
`X` (the m/z variables), using `control` as the partial variable (the per-acquisition
TIC by default). For variables ``X`` and ``Y`` with control ``I``,

```math
r_{XY \\cdot I} = \\frac{r_{XY} - r_{XI}\\,r_{YI}}{\\sqrt{(1 - r_{XI}^2)(1 - r_{YI}^2)}}
```

Positive entries indicate ions that vary together (typically the same precursor),
negative entries anticorrelated ions (different precursors). Conditioning on the
TIC removes the common intensity fluctuation shared by every ion.
"""
function partial_correlation(X::AbstractMatrix{<:Real};
                             control::AbstractVector{<:Real} = vec(sum(X, dims = 2)))
    size(X, 1) == length(control) ||
        error("partial_correlation: `control` length ($(length(control))) must equal the " *
              "number of acquisitions (rows of X = $(size(X, 1))).")
    M = size(X, 2)
    R  = cor(X)
    rc = [cor(view(X, :, j), control) for j in 1:M]
    P  = Matrix{Float64}(undef, M, M)
    @inbounds for a in 1:M, b in 1:M
        denom = sqrt((1 - rc[a]^2) * (1 - rc[b]^2))
        P[a, b] = denom == 0 ? NaN : (R[a, b] - rc[a] * rc[b]) / denom
    end
    return P
end


"""
    cmi_matrix(X; condition = vec(sum(X, dims = 2)), method = "inv", k = 3)

Symmetric `M × M` matrix of conditional mutual information ``I(X_i; X_j \\mid Z)``
between the columns of `X`, conditioned on `condition` (the per-acquisition TIC by
default). Unlike [`partial_correlation`](@ref), CMI captures nonlinear dependence
and is always ≥ 0 (it is zero only when two ions are conditionally independent).

Requires the `EntropyInvariant` package — call `using EntropyInvariant` to enable
this method (it is provided by a package extension).
"""
function cmi_matrix end
cmi_matrix(args...; kwargs...) =
    error("cmi_matrix requires the EntropyInvariant package. Run `using EntropyInvariant` " *
          "to enable it (it is provided as a package extension).")


"""
    cluster_ions(score; kind = :correlation, nclusters = nothing, linkage = :ward)

Hierarchical agglomerative clustering of the m/z variables from a square score
matrix (`partial_correlation` or `cmi_matrix`). Returns a length-`M` vector of
integer cluster labels (one per column of the abundance matrix).

`kind` selects how the score is turned into a dissimilarity: `:correlation`
(``d = 1 - r``, so positively-correlated ions cluster together) or `:cmi`
(``d = \\max(C) - C``, so high-information pairs cluster together). When
`nclusters` is `nothing` the cut level is chosen automatically by the elbow of
the within-cluster sum of squares.

Requires the `Clustering` package — call `using Clustering` to enable this method
(it is provided by a package extension).
"""
function cluster_ions end
cluster_ions(args...; kwargs...) =
    error("cluster_ions requires the Clustering package. Run `using Clustering` to enable " *
          "it (it is provided as a package extension).")


"""
    cluster_spectra(mz, X, labels)

Reconstruct one mass spectrum per cluster. `mz` are the bin centres, `X` the
`N × M` abundance matrix from [`abundance_matrix`](@ref), and `labels` the
per-column cluster assignment from [`cluster_ions`](@ref). Each output
[`MSscans`](@ref) holds the bins of one cluster with intensities set to the mean
abundance of each bin across the acquisitions, sorted by m/z. Returns the spectra
ordered by cluster label.
"""
function cluster_spectra(mz::AbstractVector{<:Real}, X::AbstractMatrix{<:Real},
                         labels::AbstractVector{<:Integer})
    (length(mz) == size(X, 2) == length(labels)) ||
        error("cluster_spectra: length(mz)=$(length(mz)), columns(X)=$(size(X, 2)) and " *
              "length(labels)=$(length(labels)) must all be equal.")
    out = MSscans[]
    for c in sort(unique(labels))
        cols  = findall(==(c), labels)
        cmz   = collect(Float64, mz[cols])
        cint  = vec(mean(view(X, :, cols), dims = 1))
        order = sortperm(cmz)
        cmz, cint = cmz[order], cint[order]
        bpi, bpidx = isempty(cint) ? (0.0, 0) : findmax(cint)
        bpm  = bpidx == 0 ? 0.0 : cmz[bpidx]
        push!(out, MSscans(Int(c), 0.0, sum(cint), cmz, cint, 2,
                           bpm, bpi, 0.0, "", "", 0.0))
    end
    return out
end
