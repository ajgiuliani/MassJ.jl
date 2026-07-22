"""
Association maps and the profile-resolution *pairwise* route.

Where [`correspondence`](@ref) (in `chimeric.jl`) *groups* co-varying ions into
per-precursor clusters, this file provides the complementary route: pick individual
correlated fragment-ion **pairs** off a two-dimensional association map, integrate
each pair's covariance **volume**, and rate it by a leave-one-out jackknife
**signal-to-noise**. This is the covariance-mapping style of analysis, done on simple
(optionally TIC-normalised) covariance so it carries no conditioning model.
"""

export MapFeature, cut_autocorrelation, map_features, covariance_features


"""
    cut_autocorrelation(A; halfwidth = 1, value = 0.0)

Return a copy of the square association map `A` with the autocorrelation diagonal
band (`|i − j| ≤ halfwidth`) set to `value` (0 by default). The diagonal of a
covariance map holds each ion's variance and dwarfs the off-diagonal cross-peaks;
suppressing that band lets pair-picking see the real correlations.
"""
function cut_autocorrelation(A::AbstractMatrix{<:Real}; halfwidth::Integer = 1,
                             value::Real = 0.0)
    n = size(A, 1)
    size(A, 2) == n || error("cut_autocorrelation: A must be square.")
    B = copy(A)
    @inbounds for i in 1:n, j in max(1, i - halfwidth):min(n, i + halfwidth)
        B[i, j] = value
    end
    return B
end


"""
    MapFeature

A correlated fragment-ion pair picked from a 2-D association map: the two m/z
(`mz1`, `mz2`), their map indices (`i`, `j`), the map `value` at the peak, the
windowed covariance `volume`, and the jackknife `snr` (`volume` / `snr` are `NaN`
until quantified by [`covariance_features`](@ref)).
"""
struct MapFeature
    mz1::Float64
    mz2::Float64
    i::Int
    j::Int
    value::Float64
    volume::Float64
    snr::Float64
end

function Base.show(io::IO, f::MapFeature)
    print(io, "MapFeature(", round(f.mz1, digits = 3), " ↔ ", round(f.mz2, digits = 3),
          ", value = ", round(f.value, sigdigits = 4))
    isnan(f.volume) || print(io, ", volume = ", round(f.volume, sigdigits = 4))
    isnan(f.snr)    || print(io, ", S/N = ", round(f.snr, digits = 2))
    print(io, ")")
end


# Clear a disk of radius r around (i,j) and its mirror (j,i) — sets them to -Inf so
# the greedy picker moves on to the next distinct peak.
function _clear!(W::AbstractMatrix, i::Integer, j::Integer, r::Integer, n::Integer)
    @inbounds for di in -r:r, dj in -r:r
        di * di + dj * dj <= r * r || continue
        a, b = i + di, j + dj
        (1 <= a <= n && 1 <= b <= n) && (W[a, b] = -Inf)
        a2, b2 = j + di, i + dj
        (1 <= a2 <= n && 1 <= b2 <= n) && (W[a2, b2] = -Inf)
    end
    return W
end


"""
    map_features(A, mz; nfeatures = 50, clearradius = 3, minvalue = 0.0,
                 edge = 0, halfwidth = 1, cut = true)

Pick up to `nfeatures` correlated pairs from the square association map `A` (its
columns indexed by `mz`) by greedy 2-D peak detection: repeatedly take the largest
remaining upper-triangle value, record it as a [`MapFeature`](@ref), and clear a disk
of radius `clearradius` around it (and its mirror). Stops at `nfeatures` or once the
largest remaining value falls to `minvalue`. With `cut = true` the autocorrelation
band (half-width `halfwidth`) is suppressed first; `edge` skips peaks within that many
pixels of the border (so an integration window fits later). Works on any symmetric
map (covariance, `partial_correlation`, `cmi_matrix`). Returns features by descending
value with `volume`/`snr` left `NaN` — fill them via [`covariance_features`](@ref).
"""
function map_features(A::AbstractMatrix{<:Real}, mz::AbstractVector{<:Real};
                      nfeatures::Integer = 50, clearradius::Integer = 3,
                      minvalue::Real = 0.0, edge::Integer = 0,
                      halfwidth::Integer = 1, cut::Bool = true)
    n = size(A, 1)
    size(A, 2) == n || error("map_features: A must be square.")
    length(mz) == n || error("map_features: length(mz) must equal size(A, 1).")
    W = cut ? cut_autocorrelation(A; halfwidth = halfwidth) : Matrix{Float64}(A)
    @inbounds for i in 1:n, j in 1:i          # keep the strict upper triangle only
        W[i, j] = -Inf
    end
    feats = MapFeature[]
    while length(feats) < nfeatures
        v, idx = findmax(W)
        (isfinite(v) && v > minvalue) || break
        i, j = Tuple(idx)
        if edge > 0 && (i <= edge || j <= edge || i > n - edge || j > n - edge)
            _clear!(W, i, j, clearradius, n)   # too close to the border: skip, keep going
            continue
        end
        push!(feats, MapFeature(mz[i], mz[j], i, j, A[i, j], NaN, NaN))
        _clear!(W, i, j, clearradius, n)
    end
    return feats
end


# Windowed covariance between column-sets from sufficient statistics.
_covwin(G, sr, sc, n) = (G .- sr * sc' ./ n) ./ (n - 1)
# Volume = sum of positive covariance over the window (a fixed-window 2-D integral).
_winvolume(cw) = sum(x -> max(x, 0.0), cw)

# Volume + jackknife S/N of the covariance peak at (i,j) over a ±window box.
function _quantify(X::AbstractMatrix{<:Real}, i::Integer, j::Integer,
                   w::Integer, jackknife::Bool)
    N, M = size(X)
    R = clamp(i - w, 1, M):clamp(i + w, 1, M)
    C = clamp(j - w, 1, M):clamp(j + w, 1, M)
    XR = X[:, R]; XC = X[:, C]
    sR = vec(sum(XR, dims = 1)); sC = vec(sum(XC, dims = 1))
    G  = XR' * XC
    vol = _winvolume(_covwin(G, sR, sC, N))
    jackknife || return (vol, NaN)
    acc = 0.0; acc2 = 0.0
    @inbounds for k in 1:N
        xr = XR[k, :]; xc = XC[k, :]
        vk = _winvolume(_covwin(G .- xr * xc', sR .- xr, sC .- xc, N - 1))
        acc += vk; acc2 += vk * vk
    end
    v̄  = acc / N
    se = sqrt(max((N - 1) / N * (acc2 - N * v̄^2), 0.0))
    return (vol, se == 0 ? (vol == 0 ? 0.0 : Inf) : vol / se)
end


"""
    covariance_features(spectra; nfeatures = 50, binsize = 1.0, threshold = 0.01,
                        mzrange = nothing, normalize = :none, window = 3,
                        clearradius = 3, halfwidth = 1, edge = window,
                        minvalue = 0.0, jackknife = true)

Profile-resolution pairwise route, end to end: build the covariance map of a chimeric
spectral series, pick correlated fragment pairs ([`map_features`](@ref)), and quantify
each by its windowed 2-D covariance `volume` and a leave-one-out jackknife `snr`.

`normalize = :tic` divides each spectrum by its total ion current before the
covariance — a simple, unconditioned way to damp a common intensity fluctuation;
`:none` uses raw intensities. `window` is the half-width (pixels) of the integration
and jackknife box. Returns a `Vector{MapFeature}` sorted by `snr` (or by `value` when
`jackknife = false`). For statistically-gated *grouping* rather than pairs, use
[`correspondence`](@ref).
"""
function covariance_features(spectra::AbstractVector{MSscans};
                             nfeatures::Integer = 50, binsize::Real = 1.0,
                             threshold::Real = 0.01, mzrange = nothing,
                             normalize::Symbol = :none, window::Integer = 3,
                             clearradius::Integer = 3, halfwidth::Integer = 1,
                             edge::Integer = -1, minvalue::Real = 0.0,
                             jackknife::Bool = true)
    am = abundance_matrix(spectra; binsize = binsize, threshold = threshold, mzrange = mzrange)
    X, mz = am.matrix, am.mz
    if normalize === :tic
        X = X ./ ifelse.(am.tic .== 0, 1.0, am.tic)
    elseif normalize !== :none
        error("covariance_features: normalize must be :none or :tic.")
    end
    A = covariance_matrix(X)
    e = edge < 0 ? window : edge
    feats = map_features(A, mz; nfeatures = nfeatures, clearradius = clearradius,
                         minvalue = minvalue, edge = e, halfwidth = halfwidth, cut = true)
    out = Vector{MapFeature}(undef, length(feats))
    for (k, f) in enumerate(feats)
        vol, snr = _quantify(X, f.i, f.j, window, jackknife)
        out[k] = MapFeature(f.mz1, f.mz2, f.i, f.j, f.value, vol, snr)
    end
    sort!(out; by = f -> jackknife ? -f.snr : -f.value)
    return out
end
