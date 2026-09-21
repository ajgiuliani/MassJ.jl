module MassJEntropyInvariantExt

using MassJ
using EntropyInvariant: conditional_mutual_information, mutual_information
import ProgressMeter

# Each variable is a single 1-D signal of N samples, reshaped to an N×1 matrix as
# EntropyInvariant's estimators expect (samples along the rows).
_ascol(v::AbstractVector{<:Real}) = reshape(collect(Float64, v), length(v), 1)

# EntropyInvariant's invariant-measure step drops zeros then requires ≥ 2 values, so a
# column populated in fewer than two acquisitions has no estimable information and would
# throw. `abundance_matrix(...; minobs)` removes such columns up front; this guards direct
# calls on an arbitrary matrix.
_estimable(col) = count(!iszero, col) >= 2

# Fill a symmetric M×M score matrix by applying `estimate(a, b)` to every upper-triangle
# ion pair. Pairs touching a non-estimable column short-circuit to 0 (MI/CMI of a
# near-constant variable is 0) rather than calling the estimator; a non-finite estimate
# (e.g. two identical non-zero values → zero spacing) is likewise mapped to 0. The pairs
# are independent (each writes its own C[a,b]/C[b,a]), so with `parallel = true` and more
# than one Julia thread they run across `Threads.@threads`; each `estimate` call stays
# serial and must be reentrant for the threaded path to be race-free.
function _pairwise_matrix(M::Int, estimate, valid::AbstractVector{Bool};
                          desc::String, progress::Bool, parallel::Bool)
    C = zeros(Float64, M, M)
    pairs = Tuple{Int,Int}[(a, b) for a in 1:M for b in (a + 1):M]
    p = progress ? ProgressMeter.Progress(length(pairs); desc = desc) : nothing
    score(a, b) = (valid[a] && valid[b]) ? (v = estimate(a, b); isfinite(v) ? v : 0.0) : 0.0
    if parallel && Threads.nthreads() > 1
        lk = ReentrantLock()                       # serialize only the progress tick
        Threads.@threads for idx in eachindex(pairs)
            a, b = pairs[idx]
            v = score(a, b)
            C[a, b] = v
            C[b, a] = v
            p === nothing || lock(() -> ProgressMeter.next!(p), lk)
        end
    else
        for (a, b) in pairs
            v = score(a, b)
            C[a, b] = v
            C[b, a] = v
            p === nothing || ProgressMeter.next!(p)
        end
    end
    return C
end

function MassJ.mi(X::AbstractMatrix{<:Real}; method::String = "inv", k::Int = 3,
                  progress::Bool = true, parallel::Bool = true)
    M     = size(X, 2)
    cols  = [_ascol(view(X, :, j)) for j in 1:M]
    valid = [_estimable(view(X, :, j)) for j in 1:M]
    return _pairwise_matrix(M,
        (a, b) -> mutual_information(cols[a], cols[b]; method = method, k = k), valid;
        desc = "mi: ", progress = progress, parallel = parallel)
end

function MassJ.cmi(X::AbstractMatrix{<:Real}, control::AbstractVecOrMat{<:Real};
                   method::String = "inv", k::Int = 3,
                   progress::Bool = true, parallel::Bool = true)
    size(X, 1) == size(control, 1) ||
        error("cmi: control must have $(size(X, 1)) rows (one per acquisition), " *
              "got $(size(control, 1)).")
    control isa AbstractMatrix && size(control, 2) > 1 &&
        error("cmi conditions on a single variable — EntropyInvariant's conditional MI takes a " *
              "1-D conditioning variable, so pass a vector. (For multivariate conditioning use pcor.)")
    zc = control isa AbstractMatrix ? vec(control) : control
    _estimable(zc) ||
        error("cmi: control must have at least 2 non-zero values to condition on.")
    M     = size(X, 2)
    Z     = _ascol(zc)
    cols  = [_ascol(view(X, :, j)) for j in 1:M]
    valid = [_estimable(view(X, :, j)) for j in 1:M]
    return _pairwise_matrix(M,
        (a, b) -> conditional_mutual_information(cols[a], cols[b], Z; method = method, k = k), valid;
        desc = "cmi: ", progress = progress, parallel = parallel)
end

end # module
