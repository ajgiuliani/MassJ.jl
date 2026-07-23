"""
Statistical decomposition of chimeric (multiplexed) tandem mass spectra.

Implements the data-driven workflow of Truong, Nahon & Giuliani
(*J. Am. Soc. Mass Spectrom.* **2026**, 37, 649): a series of `N` tandem mass
spectra recorded during photon activation is reduced to an abundance matrix
(rows = acquisitions, columns = binned m/z variables); the columns are scored
against one another by partial Pearson correlation ([`pcor`](@ref))
or conditional mutual information ([`cmi`](@ref)); hierarchical clustering
of the score matrix ([`cluster_ions`](@ref)) groups the ions by precursor, and
[`cluster_spectra`](@ref) reconstructs one mass spectrum per cluster.

[`cmi`](@ref) requires the `EntropyInvariant` package and
[`cluster_ions`](@ref) the `Clustering` package; both are loaded on demand as
package extensions (`using EntropyInvariant` / `using Clustering`).
"""

export abundance_matrix, cluster_ions, cluster_spectra


"""
    abundance_matrix(spectra; binsize = 1.0, threshold = 0.01, mzrange = nothing,
                     drop_constant = true, minobs = 2)

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

When `drop_constant = true` (default) bins that carry no usable signal are removed:
those with zero variance across acquisitions, and those populated in fewer than
`minobs` acquisitions. A fragment seen in a single scan carries no co-variation for
any of the scores, and fewer than two non-zero observations also breaks the
information-theoretic estimators ([`mi`](@ref)/[`cmi`](@ref)); raise `minobs` to also
drop rarely-seen bins.
"""
function abundance_matrix(spectra::AbstractVector{MSscans};
                          binsize::Real = 1.0,
                          threshold::Real = 0.01,
                          mzrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
                          drop_constant::Bool = true,
                          minobs::Integer = 2)
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
        keep = [var(@view X[:, j]) > 0 && count(!iszero, @view X[:, j]) >= minobs
                for j in 1:M]
        X = X[:, keep]
        centers = centers[keep]
    end

    tic = vec(sum(X, dims = 2))
    return (matrix = X, mz = centers, tic = tic)
end


"""
    cluster_ions(score; kind = :correlation, nclusters = nothing, linkage = :ward)

Hierarchical agglomerative clustering of the m/z variables from a square score
matrix (`pcor` or `cmi`). Returns a length-`M` vector of
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
#  Score functions (cov / cor / pcor / mi / cmi), leave-one-out and permutation
#  significance, FDR gating, and the hands-off correspondence pipeline. `cov`/`cor`
#  reuse Statistics; `pcor` is covariance-based; `mi`/`cmi` come from EntropyInvariant
#  (package extension). All are pure: matrix in, score matrix out — binning lives in
#  `abundance_matrix`, significance in `jackknife`/`permutation`.
# =====================================================================

export cov, cor, pcor, mi, cmi, jackknife, permutation, fdr_adjust, correspondence

# --- cov / cor: reuse Statistics on a matrix; add an MSscans convenience ---
# (recommended path: `abundance_matrix` first, then pass the matrix; the MSscans form
#  runs default `abundance_matrix` and is a quick, less-recommended shortcut.)
Statistics.cov(s::AbstractVector{MSscans}) = cov(abundance_matrix(s).matrix)
Statistics.cor(s::AbstractVector{MSscans}) = cor(abundance_matrix(s).matrix)

"""
    pcor(X, control)                 # X: matrix (rows = acquisitions)
    pcor(spectra[, control])         # MSscans series (bins with abundance_matrix first)

Symmetric `M × M` partial (Pearson) correlation between the columns of `X`, conditioned
on `control` — a length-`N` vector (one conditioning variable) or an `N × p` matrix
(several). Computed from the correlations of the residuals after regressing every column
on `control` (with an intercept), so it generalises cleanly to multivariate control.
Pass `control = am.tic` to condition on the total ion current; `cor(X)` is the
unconditioned Pearson correlation.
"""
function pcor(X::AbstractMatrix{<:Real}, control::AbstractVecOrMat{<:Real})
    N = size(X, 1)
    size(control, 1) == N ||
        error("pcor: control must have $N rows (one per acquisition), got $(size(control, 1)).")
    Xf = float.(X)
    Z  = hcat(ones(N), control)              # design matrix with intercept
    R  = Xf .- Z * (Z \ Xf)                  # residuals after regressing out control
    C  = cor(R)
    C[.!isfinite.(C)] .= 0.0                 # constant residual column -> 0
    return C
end
pcor(s::AbstractVector{MSscans}) = (am = abundance_matrix(s); pcor(am.matrix, am.tic))
pcor(s::AbstractVector{MSscans}, control) = pcor(abundance_matrix(s).matrix, control)

"""
    mi(X; method = "inv", k = 3, progress = true, parallel = true)

Symmetric `M × M` (unconditioned) mutual information between the columns of `X` — the
information-theoretic sibling of `cor`: nonlinear and always ≥ 0. Requires the
`EntropyInvariant` package (package extension); [`cmi`](@ref) is its conditioned form.

The `M(M-1)/2` pair estimates are independent, so with `parallel = true` (the default)
they run across the threads Julia was started with (`julia -t N`); each estimate stays
serial. Pass `parallel = false` to keep the pair loop single-threaded — e.g. when you
parallelize an enclosing [`jackknife`](@ref)/[`permutation`](@ref) loop instead, to avoid
nesting the two.
"""
function mi end
mi(args...; kwargs...) =
    error("mi requires the EntropyInvariant package. Run `using EntropyInvariant` to enable it.")
mi(s::AbstractVector{MSscans}; kwargs...) = mi(abundance_matrix(s).matrix; kwargs...)

"""
    cmi(X, control; method = "inv", k = 3, progress = true, parallel = true)

Symmetric `M × M` conditional mutual information ``I(X_i; X_j ∣ Z)`` between the columns
of `X`, conditioned on `control` (a length-`N` vector — a single conditioning variable;
EntropyInvariant's conditional MI is 1-D). The nonlinear counterpart of [`pcor`](@ref),
always ≥ 0 (up to finite-sample estimator noise). Requires the `EntropyInvariant` package.

Like [`mi`](@ref), the pair estimates run across Julia's threads when `parallel = true`
(the default); pass `parallel = false` to keep the loop serial.
"""
function cmi end
cmi(args...; kwargs...) =
    error("cmi requires the EntropyInvariant package. Run `using EntropyInvariant` to enable it.")
cmi(s::AbstractVector{MSscans}; kwargs...) =
    (am = abundance_matrix(s); cmi(am.matrix, am.tic; kwargs...))
cmi(s::AbstractVector{MSscans}, control; kwargs...) =
    cmi(abundance_matrix(s).matrix, control; kwargs...)

# Apply a score to (X[, control]); the `q` variants silence a slow score's own progress
# bar (mi/cmi) when they are re-run inside a resampling loop.
_apply(f,  X, ::Nothing) = f(X)
_apply(f,  X, control)   = f(X, control)
_applyq(f, X, ::Nothing) = (f === mi || f === cmi) ? f(X; progress = false) : f(X)
_applyq(f, X, control)   = (f === mi || f === cmi) ? f(X, control; progress = false) : f(X, control)
_subrows(c::AbstractVector, mask) = view(c, mask)
_subrows(c::AbstractMatrix, mask) = view(c, mask, :)


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


"""
    jackknife(score, X[, control]) -> (estimate, pvalue)

Leave-one-out jackknife significance of a `score` (`cov`, `cor`, or `pcor`) between the
columns of the abundance matrix `X`. Recomputes `score` with each acquisition removed
and turns the spread into a two-sided p-value via the z-score `estimate / se`
(`se = √((N-1)/N·Σ(θₖ−θ̄)²)`); a small p marks a pair whose association is stable rather
than driven by a single outlier acquisition. Pass `control` for a conditioned score
(`pcor`) — it is subset alongside `X`. Not valid for `mi`/`cmi` (their null is not
zero-centred) — use [`permutation`](@ref). `progress` shows a bar over the `N` recomputes.
"""
function jackknife(score, X::AbstractMatrix{<:Real}, control = nothing; progress::Bool = true)
    (score === mi || score === cmi) &&
        error("jackknife is not valid for mi/cmi (their null is not zero-centred); use permutation.")
    N, M = size(X)
    N > 2 || error("jackknife: need at least three acquisitions.")
    est  = _apply(score, X, control)
    acc  = zeros(Float64, M, M); acc2 = zeros(Float64, M, M)
    mask = trues(N)
    p = progress ? Progress(N; desc = "jackknife: ") : nothing
    for k in 1:N
        mask[k] = false
        θ = _apply(score, view(X, mask, :), control === nothing ? nothing : _subrows(control, mask))
        acc .+= θ; acc2 .+= θ .^ 2
        mask[k] = true
        p === nothing || ProgressMeter.next!(p)
    end
    θbar = acc ./ N
    se   = sqrt.(max.(((N - 1) / N) .* (acc2 .- N .* θbar .^ 2), 0.0))
    return (estimate = est, pvalue = _z_to_p.(est ./ se))
end

"""
    permutation(score, X[, control]; nperm = 200, rng = Random.default_rng()) -> (estimate, pvalue)

Permutation significance of a `score` (any of `cov`/`cor`/`pcor`/`mi`/`cmi`). Builds a
null by independently permuting each column of `X` (holding `control` fixed), recomputes
`score` `nperm` times, and reports `p = (1 + #{null ⋈ obs}) / (nperm + 1)`. This is the
route for `mi`/`cmi`, whose null is not zero-centred. `progress` shows a bar over `nperm`.
"""
function permutation(score, X::AbstractMatrix{<:Real}, control = nothing;
                     nperm::Integer = 200, rng = Random.default_rng(), progress::Bool = true)
    N, M   = size(X)
    signed = !(score === mi || score === cmi)     # signed statistics compare on |·|
    est    = _applyq(score, X, control)
    cmpobs = signed ? abs.(est) : est
    cnt = zeros(Int, M, M)
    Xp  = Matrix{Float64}(undef, N, M)
    p = progress ? Progress(nperm; desc = "permutation: ") : nothing
    for _ in 1:nperm
        for j in 1:M
            col = @view Xp[:, j]; copyto!(col, @view X[:, j]); shuffle!(rng, col)
        end
        ns   = _applyq(score, Xp, control)
        cmpn = signed ? abs.(ns) : ns
        cnt .+= (cmpn .>= cmpobs)
        p === nothing || ProgressMeter.next!(p)
    end
    return (estimate = est, pvalue = (1 .+ cnt) ./ (nperm + 1))
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
    correspondence(spectra; measure = :pcor, control = :tic, method = :jackknife, q = 0.05,
                   binsize = 1.0, threshold = 0.01, mzrange = nothing, nperm = 200,
                   rng = Random.default_rng(), nclusters = nothing, linkage = :ward,
                   progress = true)

Hands-off decomposition of a chimeric / multiplexed spectral series into per-precursor
spectra — a thin chain over the explicit steps: [`abundance_matrix`](@ref) → the
`measure` score (`:cov`, `:cor`, `:pcor`, `:mi`, `:cmi`) → significance
(`method = :jackknife` or `:permutation`) → Benjamini-Hochberg FDR gate at level `q`
([`fdr_adjust`](@ref)) → [`cluster_ions`](@ref) → [`cluster_spectra`](@ref).

`control` sets the conditioning variable for `:pcor`/`:cmi`: `:tic` (default, the
per-acquisition TIC), `:none`, or a per-acquisition vector. `:mi`/`:cmi` require
`method = :permutation` and `EntropyInvariant`; clustering requires `Clustering`. If no
edge survives the gate, the ungated score is used (with a warning). Returns
`(labels, spectra, score, pvalue, significant, mz)`.
"""
function correspondence(spectra::AbstractVector{MSscans};
                        measure::Symbol = :pcor, control = :tic, method::Symbol = :jackknife,
                        q::Real = 0.05, binsize::Real = 1.0, threshold::Real = 0.01,
                        mzrange = nothing, nperm::Integer = 200, rng = Random.default_rng(),
                        nclusters = nothing, linkage::Symbol = :ward, progress::Bool = true)
    scorefn = measure === :cov  ? cov : measure === :cor ? cor : measure === :pcor ? pcor :
              measure === :mi   ? mi  : measure === :cmi ? cmi :
              error("correspondence: measure must be :cov, :cor, :pcor, :mi, or :cmi.")
    am    = abundance_matrix(spectra; binsize = binsize, threshold = threshold, mzrange = mzrange)
    X, mz = am.matrix, am.mz
    conditioned = measure === :pcor || measure === :cmi
    ctrl = !conditioned ? nothing :
           control === :tic ? am.tic : control === :none ? nothing : control
    sig  = method === :jackknife  ? jackknife(scorefn, X, ctrl; progress = progress) :
           method === :permutation ? permutation(scorefn, X, ctrl; nperm = nperm, rng = rng, progress = progress) :
           error("correspondence: method must be :jackknife or :permutation.")
    S     = sig.estimate
    fg    = fdr_adjust(sig.pvalue; q = q)
    kind  = (measure === :mi || measure === :cmi) ? :cmi : :correlation
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
            score = S, pvalue = sig.pvalue, significant = fg.significant, mz = mz)
end
