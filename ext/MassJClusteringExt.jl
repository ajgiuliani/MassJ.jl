module MassJClusteringExt

using MassJ
using Clustering: hclust, cutree

# Turn a similarity-like score matrix into a symmetric, non-negative dissimilarity.
function _score_to_distance(S::AbstractMatrix{<:Real}, kind::Symbol)
    M = size(S, 1)
    D = Matrix{Float64}(undef, M, M)
    if kind === :correlation
        @inbounds for i in 1:M, j in 1:M
            D[i, j] = i == j ? 0.0 : clamp(1 - S[i, j], 0.0, 2.0)
        end
    elseif kind === :cmi
        smax = maximum(S)
        @inbounds for i in 1:M, j in 1:M
            D[i, j] = i == j ? 0.0 : max(smax - S[i, j], 0.0)
        end
    else
        error("cluster_ions: `kind` must be :correlation or :cmi (got $(kind)).")
    end
    return (D .+ D') ./ 2     # enforce exact symmetry
end

# Within-cluster sum of squared distances for a k-cut of the dendrogram.
function _wss(D::AbstractMatrix{<:Real}, h, k::Int)
    labels = cutree(h; k = k)
    s = 0.0
    @inbounds for c in unique(labels)
        idx = findall(==(c), labels)
        for a in 1:length(idx), b in (a + 1):length(idx)
            s += D[idx[a], idx[b]]^2
        end
    end
    return s
end

# Elbow ("knee") of the WSS-vs-k curve: the k whose point lies farthest from the
# straight line joining the first and last points (Kneedle-style).
function _elbow_k(D::AbstractMatrix{<:Real}, h, kmax::Int)
    kmax = max(kmax, 1)
    wss = [_wss(D, h, k) for k in 1:kmax]
    n = length(wss)
    n <= 2 && return n
    x1, y1 = 1.0, wss[1]
    x2, y2 = float(n), wss[end]
    denom = sqrt((y2 - y1)^2 + (x2 - x1)^2)
    denom == 0 && return 1
    best, bestk = -Inf, 1
    for i in 1:n
        d = abs((y2 - y1) * i - (x2 - x1) * wss[i] + x2 * y1 - y2 * x1) / denom
        if d > best
            best, bestk = d, i
        end
    end
    return bestk
end

function MassJ.cluster_ions(score::AbstractMatrix{<:Real};
                            kind::Symbol = :correlation,
                            nclusters::Union{Nothing,Integer} = nothing,
                            linkage::Symbol = :ward,
                            kmax::Integer = min(10, size(score, 1) - 1))
    size(score, 1) == size(score, 2) ||
        error("cluster_ions: `score` must be a square matrix (got $(size(score))).")
    D = _score_to_distance(score, kind)
    h = hclust(D; linkage = linkage)
    k = nclusters === nothing ? _elbow_k(D, h, Int(kmax)) : Int(nclusters)
    return cutree(h; k = k)
end

end # module
