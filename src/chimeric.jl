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


# =====================================================================
#  Conditioned-association significance, FDR gating, and the hands-off
#  correspondence pipeline. These add a per-edge statistical confidence to
#  the score matrices above (partial_correlation / cmi_matrix), so the ion
#  grouping can be gated by significance rather than a hand-picked threshold.
# =====================================================================

export covariance_matrix, jackknife_significance, permutation_significance,
       fdr_adjust, correspondence


"""
    covariance_matrix(X; corrected = true)

Symmetric `M × M` sample covariance between the columns (m/z variables) of the
abundance matrix `X` (rows = acquisitions) — the classical covariance map. The
diagonal holds each variable's variance; `corrected = false` divides by `N`
instead of `N - 1`. Unlike [`partial_correlation`](@ref) it applies no
conditioning, so a common intensity fluctuation shared by every ion appears as a
positive background; condition it out with `partial_correlation` (linear) or
`cmi_matrix` (information-theoretic) when that matters.
"""
function covariance_matrix(X::AbstractMatrix{<:Real}; corrected::Bool = true)
    N = size(X, 1)
    N > 1 || error("covariance_matrix: need at least two acquisitions.")
    Xc = X .- mean(X, dims = 1)
    return (Xc' * Xc) ./ (corrected ? (N - 1) : N)
end


# Two-sided standard-normal tail p = erfc(|z|/√2), via the Abramowitz & Stegun
# 7.1.26 erf approximation (|error| < 1.5e-7) — no SpecialFunctions dependency.
function _normal_sf2(z::Real)
    x = abs(z) / sqrt(2)
    t = 1 / (1 + 0.3275911 * x)
    erf = 1 - (((((1.061405429t - 1.453152027)t + 1.421413741)t -
                 0.284496736)t + 0.254829592)t) * exp(-x * x)
    return clamp(1 - erf, 0.0, 1.0)          # erfc(x) for x ≥ 0
end
_z_to_p(z::Real) = isnan(z) ? 1.0 : (isinf(z) ? 0.0 : _normal_sf2(z))


# M×M pairwise statistic among the first `M` columns from the sufficient
# statistics (column sums `s1`, Gram matrix `S2`) computed on `n` rows.
# For :pcorr the control occupies the last column (index length(s1)).
function _pairwise_stat(s1::AbstractVector, S2::AbstractMatrix, n::Integer,
                        statistic::Symbol, M::Integer)
    Cov = (S2 .- (s1 * s1') ./ n) ./ (n - 1)          # P×P
    statistic === :covariance && return Cov[1:M, 1:M]
    sd = sqrt.(max.(diag(Cov), 0.0))
    R = Cov ./ (sd * sd')
    R[.!isfinite.(R)] .= 0.0
    statistic === :correlation && return R[1:M, 1:M]
    c  = length(s1)                                    # control column
    rc = R[1:M, c]
    denom = sqrt.((1 .- rc .^ 2) * (1 .- rc .^ 2)')
    P = (R[1:M, 1:M] .- rc * rc') ./ denom
    P[.!isfinite.(P)] .= 0.0
    return P
end


"""
    jackknife_significance(X; statistic = :pcorr, control = vec(sum(X, dims = 2)))

Delete-one jackknife significance of a pairwise association between the columns of
the abundance matrix `X`. Returns a named tuple `(estimate, se, z)` of `M × M`
matrices: the full-sample statistic, its jackknife standard error
(`√((N-1)/N · Σ(θₖ − θ̄)²)` over the leave-one-out replicates `θₖ`), and the
z-score `estimate / se`. A large `|z|` marks a pair whose association is stable
under leave-one-out resampling — i.e. supported by the whole series rather than a
single outlier acquisition, the failure mode that inflates raw covariance.

`statistic` is `:covariance`, `:correlation`, or `:pcorr` (partial correlation
conditioned on `control`, the per-acquisition TIC by default). Convert the z-score
to a two-sided p-value with the standard-normal tail. Cost is `O(N·M²)`; for very
large `M` (profile-resolution maps) prefer the windowed feature route.
"""
function jackknife_significance(X::AbstractMatrix{<:Real};
                                statistic::Symbol = :pcorr,
                                control::AbstractVector{<:Real} = vec(sum(X, dims = 2)))
    statistic in (:covariance, :correlation, :pcorr) ||
        error("jackknife_significance: statistic must be :covariance, :correlation, or :pcorr.")
    N, M = size(X)
    N > 2 || error("jackknife_significance: need at least three acquisitions.")
    D  = statistic === :pcorr ? hcat(float.(X), float.(control)) : float.(X)
    s1 = vec(sum(D, dims = 1))
    S2 = D' * D
    est  = _pairwise_stat(s1, S2, N, statistic, M)
    acc  = zeros(Float64, M, M)
    acc2 = zeros(Float64, M, M)
    for k in 1:N
        d   = D[k, :]
        θ   = _pairwise_stat(s1 .- d, S2 .- d * d', N - 1, statistic, M)
        acc  .+= θ
        acc2 .+= θ .^ 2
    end
    θbar = acc ./ N
    v    = max.(((N - 1) / N) .* (acc2 .- N .* θbar .^ 2), 0.0)
    se   = sqrt.(v)
    return (estimate = est, se = se, z = est ./ se)
end


"""
    permutation_significance(X; statistic = :cmi, condition = vec(sum(X, dims = 2)),
                             nperm = 200, rng = Random.default_rng(),
                             tail = :greater, kwargs...)

Permutation significance of a pairwise association between the columns of `X`.
Returns `(estimate, pvalue)` of `M × M` matrices. The null is built by
independently permuting each column of `X` (destroying inter-column association
while preserving each variable's marginal), recomputing the statistic `nperm`
times, and counting how often the null reaches the observed value:
`p = (1 + #{null ⋈ obs}) / (nperm + 1)`. This is the significance route for a
statistic whose null is not zero-centred, notably conditional mutual information
(`:cmi`, requires `EntropyInvariant`), where a jackknife z-score does not apply.
`statistic` is `:cmi`, `:pcorr`, `:correlation`, or `:covariance`; `tail` is
`:greater` (non-negative statistics like CMI) or `:two` (signed statistics).

Note: the null permutes columns unconditionally, so it tests association beyond
chance rather than a strict conditional-independence null; extra `kwargs` are
forwarded to `cmi_matrix`.
"""
function permutation_significance(X::AbstractMatrix{<:Real};
                                  statistic::Symbol = :cmi,
                                  condition::AbstractVector{<:Real} = vec(sum(X, dims = 2)),
                                  nperm::Integer = 200,
                                  rng = Random.default_rng(),
                                  tail::Symbol = :greater, kwargs...)
    N, M = size(X)
    statfun =
        statistic === :cmi         ? (Xm -> cmi_matrix(Xm; condition = condition, kwargs...)) :
        statistic === :pcorr       ? (Xm -> partial_correlation(Xm; control = condition))     :
        statistic === :correlation ? (Xm -> cor(Xm))                                          :
        statistic === :covariance  ? (Xm -> covariance_matrix(Xm))                            :
        error("permutation_significance: statistic must be :cmi, :pcorr, :correlation, or :covariance.")
    obs    = statfun(X)
    cmpobs = tail === :two ? abs.(obs) : obs
    cnt    = zeros(Int, M, M)
    Xp     = Matrix{Float64}(undef, N, M)
    for _ in 1:nperm
        for j in 1:M
            col = @view Xp[:, j]
            copyto!(col, @view X[:, j])
            shuffle!(rng, col)
        end
        ns   = statfun(Xp)
        cmpn = tail === :two ? abs.(ns) : ns
        cnt .+= (cmpn .>= cmpobs)
    end
    return (estimate = obs, pvalue = (1 .+ cnt) ./ (nperm + 1))
end


"""
    fdr_adjust(P; q = 0.05)

Benjamini-Hochberg false-discovery-rate control over the pairwise p-value matrix
`P` (symmetric, `M × M`). Only the unique upper-triangle pairs (`i < j`) are tested;
the diagonal is never significant. Returns `(significant, qvalue)`: a symmetric
`BitMatrix` of edges whose BH q-value is ≤ `q`, and the symmetric matrix of BH
q-values (1.0 on untested cells). Use it to keep only edges whose association
survives multiple-comparison correction before clustering.
"""
function fdr_adjust(P::AbstractMatrix{<:Real}; q::Real = 0.05)
    M   = size(P, 1)
    idx = [(i, j) for i in 1:M for j in i+1:M]
    pv  = [P[i, j] for (i, j) in idx]
    m   = length(pv)
    m == 0 && return (significant = falses(M, M), qvalue = fill(1.0, M, M))
    ord  = sortperm(pv)
    qval = fill(1.0, m)
    minq = 1.0
    for r in m:-1:1
        k    = ord[r]
        minq = min(minq, pv[k] * m / r)
        qval[k] = min(minq, 1.0)
    end
    sig = falses(M, M)
    Q   = fill(1.0, M, M)
    for (t, (i, j)) in enumerate(idx)
        Q[i, j] = Q[j, i] = qval[t]
        if qval[t] <= q
            sig[i, j] = sig[j, i] = true
        end
    end
    return (significant = sig, qvalue = Q)
end


"""
    correspondence(spectra; score = :pcorr, significance = :jackknife, q = 0.05,
                   binsize = 1.0, threshold = 0.01, mzrange = nothing,
                   control = nothing, nperm = 200, rng = Random.default_rng(),
                   nclusters = nothing, linkage = :ward, kwargs...)

Hands-off decomposition of a chimeric / multiplexed spectral series into
per-precursor spectra, gating the association graph by statistical significance so
the grouping needs no manually chosen correlation threshold.

Pipeline: [`abundance_matrix`](@ref) → association score (`:pcorr`, `:covariance`,
or `:cmi`) with its per-edge significance (`:jackknife` z→p for the correlation
family, `:permutation` for `:cmi`) → Benjamini-Hochberg FDR gate at level `q`
([`fdr_adjust`](@ref)) → [`cluster_ions`](@ref) on the significant-edge graph →
[`cluster_spectra`](@ref). `control` defaults to the per-acquisition TIC. Requires
the `Clustering` package (and `EntropyInvariant` when `score = :cmi`).

Returns `(labels, spectra, score, pvalue, significant, mz)`.
"""
function correspondence(spectra::AbstractVector{MSscans};
                        score::Symbol = :pcorr,
                        significance::Symbol = (score === :cmi ? :permutation : :jackknife),
                        q::Real = 0.05, binsize::Real = 1.0, threshold::Real = 0.01,
                        mzrange = nothing, control = nothing,
                        nperm::Integer = 200, rng = Random.default_rng(),
                        nclusters = nothing, linkage::Symbol = :ward, kwargs...)
    am   = abundance_matrix(spectra; binsize = binsize, threshold = threshold, mzrange = mzrange)
    X, mz = am.matrix, am.mz
    ctrl = control === nothing ? am.tic : control

    S, kind =
        score === :pcorr      ? (partial_correlation(X; control = ctrl), :correlation) :
        score === :covariance ? (covariance_matrix(X),                   :correlation) :
        score === :cmi        ? (cmi_matrix(X; condition = ctrl, kwargs...), :cmi)      :
        error("correspondence: score must be :pcorr, :covariance, or :cmi.")

    p = if significance === :jackknife
        _z_to_p.(jackknife_significance(X; statistic = score, control = ctrl).z)
    elseif significance === :permutation
        permutation_significance(X; statistic = score, condition = ctrl, nperm = nperm,
                                 rng = rng, tail = (score === :cmi ? :greater : :two),
                                 kwargs...).pvalue
    else
        error("correspondence: significance must be :jackknife or :permutation.")
    end

    fg = fdr_adjust(p; q = q)
    gated = if any(fg.significant)
        g = S .* fg.significant
        @inbounds for i in axes(g, 1)          # restore self-similarity for clustering
            g[i, i] = S[i, i]
        end
        g
    else
        @warn "correspondence: no edges survived FDR gating at q=$q; clustering on the " *
              "ungated score. Try a larger nperm, a looser q, or more acquisitions."
        copy(S)
    end
    labels = cluster_ions(gated; kind = kind, nclusters = nclusters, linkage = linkage)
    return (labels = labels, spectra = cluster_spectra(mz, X, labels),
            score = S, pvalue = p, significant = fg.significant, mz = mz)
end
