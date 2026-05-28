module MassJEntropyInvariantExt

using MassJ
using EntropyInvariant: conditional_mutual_information

# Each variable is a single 1-D signal of N samples. EntropyInvariant's matrix
# methods expect samples along the rows (`dim = 1`, one column per variable), so
# a univariate signal is reshaped to an N×1 matrix.
_ascol(v::AbstractVector{<:Real}) = reshape(collect(Float64, v), length(v), 1)

function MassJ.cmi_matrix(X::AbstractMatrix{<:Real};
                          condition::AbstractVector{<:Real} = vec(sum(X, dims = 2)),
                          method::String = "inv", k::Int = 3)
    size(X, 1) == length(condition) ||
        error("cmi_matrix: `condition` length ($(length(condition))) must equal the " *
              "number of acquisitions (rows of X = $(size(X, 1))).")
    M    = size(X, 2)
    Z    = _ascol(condition)
    cols = [_ascol(view(X, :, j)) for j in 1:M]
    C    = zeros(Float64, M, M)
    @inbounds for a in 1:M, b in (a + 1):M
        v = conditional_mutual_information(cols[a], cols[b], Z; method = method, k = k)
        C[a, b] = v
        C[b, a] = v
    end
    return C
end

end # module
