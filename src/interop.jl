# Bridges between MassJ's typed spectrum containers and the wider Julia data /
# machine-learning ecosystem. The pattern is the same as the Tables / Measurements
# / Unitful extensions: keep the typed core, expose thin converters to a feature
# *matrix* (the universal ML interface) and back. Results from any external
# algorithm (outlier scores, cluster labels, …) are one value per row, and the
# matrix is row-aligned to the input series, so they map straight back by index.

export featurize, select_spectra, from_matrix, spectra_table


"""
    featurize(scans; method = :bins, binsize = 1.0, mzrange = nothing,
              threshold = 0.0, drop_constant = false) -> (matrix, mz)

Reduce a series of spectra (`AbstractVector{MSscans}`, e.g. an [`MSrun`](@ref) or
the result of [`load`](@ref)) to a dense feature matrix suitable for clustering,
outlier detection, random forests, PCA, … `matrix` is `N×M` (one **row per
spectrum**, row-aligned to the input; one column per m/z bin) and `mz` holds the
`M` bin centres.

`method = :bins` (currently the only method) sums each spectrum's intensities into
`binsize`-wide m/z bins on a common axis spanning `mzrange` (or the global m/z
range). `threshold` (fraction of each spectrum's base peak) and `drop_constant`
(remove zero-variance bins) default to keeping *everything*, so the matrix is a
faithful binned image of the spectra that [`from_matrix`](@ref) can rebuild.

```julia
scans   = load("series/")
X, mz   = featurize(scans; binsize = 0.5)     # N×M matrix + bin centres
# … run any ecosystem algorithm on X (rows = spectra) …
```

See [`spectra_table`](@ref) for a [Tables.jl](https://github.com/JuliaData/Tables.jl)
view (e.g. `DataFrame`), [`select_spectra`](@ref) to map a per-row mask back to a
`Vector{MSscans}`, and [`from_matrix`](@ref) to rebuild spectra from a modified
matrix.
"""
function featurize(scans::AbstractVector{MSscans};
                   method::Symbol = :bins,
                   binsize::Real = 1.0,
                   mzrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
                   threshold::Real = 0.0,
                   drop_constant::Bool = false)
    method === :bins ||
        error("featurize: only method = :bins is implemented (got :$method)")
    A = abundance_matrix(scans; binsize = binsize, threshold = threshold,
                         mzrange = mzrange, drop_constant = drop_constant)
    return (matrix = A.matrix, mz = A.mz)
end


"""
    select_spectra(scans, mask::AbstractVector{Bool}) -> Vector{MSscans}
    select_spectra(scans, idx::AbstractVector{<:Integer}) -> Vector{MSscans}

Select a subset of a spectrum series by a per-row boolean `mask` (e.g. an inlier
mask from an outlier analysis run on [`featurize`](@ref)`(scans)`) or by integer
indices. Equivalent to `scans[mask]` / `scans[idx]` with a length check; named for
discoverability in the round-trip workflow. (`select_spectra`, not `select`, to
avoid clashing with `DataFrames.select`.)
"""
function select_spectra(scans::AbstractVector{MSscans}, mask::AbstractVector{Bool})
    length(mask) == length(scans) ||
        error("select_spectra: length(mask) ($(length(mask))) != length(scans) ($(length(scans)))")
    return scans[mask]
end
select_spectra(scans::AbstractVector{MSscans}, idx::AbstractVector{<:Integer}) = scans[idx]


# Scalar value of a possibly-vector provenance field, or a default when empty.
_one(v::AbstractVector, default) = isempty(v) ? default : first(v)

# Build one MSscans from an m/z axis + intensity vector, carrying a template's
# per-scan metadata when supplied.
function _spectrum_from_vectors(mz::Vector{Float64}, int::Vector{Float64},
                                t::Union{MSscans,Nothing})
    bpi = isempty(int) ? 0.0 : maximum(int)
    bpm = isempty(int) ? 0.0 : mz[argmax(int)]
    tic = isempty(int) ? 0.0 : sum(int)
    if t === nothing
        return MSscans(1, 0.0, tic, mz, int, 1, bpm, bpi, 0.0, "", "", 0.0)
    end
    return MSscans(_one(t.num, 1), _one(t.rt, 0.0), tic, mz, int, _one(t.level, 1),
                   bpm, bpi, _one(t.precursor, 0.0), _one(t.polarity, ""),
                   _one(t.activationMethod, ""), _one(t.collisionEnergy, 0.0),
                   _one(t.chargeState, 0), t.spectrumType, _one(t.driftTime, -1.0),
                   _one(t.compensationVoltage, 0.0), t.mobilityType, copy(t.metadata))
end


"""
    from_matrix(X::AbstractMatrix, mz::AbstractVector; template = nothing) -> Vector{MSscans}

Rebuild a spectrum series from a feature matrix — the inverse of
[`featurize`](@ref). Each **row** of `X` becomes one [`MSscans`](@ref) on the
**binned** m/z axis `mz` (the bin centres), with intensities taken from that row.
This is the right call after transforming the intensity matrix in a way
[`select_spectra`](@ref) cannot express (scaling, denoising, imputation, an
autoencoder reconstruction, …).

`size(X, 2)` must equal `length(mz)`. Pass `template` (a `Vector{MSscans}` with one
entry per row, e.g. the original `scans`) to carry each spectrum's metadata
(retention time, scan number, polarity, MS level, …) onto the rebuilt spectrum;
without it, neutral defaults are used. The result lives on the binned axis, so it
is a *resampled* view of the originals, not a bit-identical copy.

```julia
X, mz   = featurize(scans)
X2      = denoise(X)                              # any per-row transform
rebuilt = from_matrix(X2, mz; template = scans)   # back to Vector{MSscans}
```
"""
function from_matrix(X::AbstractMatrix{<:Real}, mz::AbstractVector{<:Real};
                     template::Union{AbstractVector{MSscans},Nothing} = nothing)
    n, m = size(X)
    m == length(mz) ||
        error("from_matrix: size(X, 2) ($m) != length(mz) ($(length(mz)))")
    template === nothing || length(template) == n ||
        error("from_matrix: length(template) ($(length(template))) != number of rows ($n)")
    axis = collect(Float64, mz)
    out = Vector{MSscans}(undef, n)
    for i in 1:n
        t = template === nothing ? nothing : template[i]
        out[i] = _spectrum_from_vectors(copy(axis), Float64.(@view X[i, :]), t)
    end
    return out
end


"""
    spectra_table(scans; binsize = 1.0, mzrange = nothing, metadata = true) -> NamedTuple

A [Tables.jl](https://github.com/JuliaData/Tables.jl)-ready view of a spectrum
series: a `NamedTuple` of equal-length vectors with one **row per spectrum**.
Columns are the per-spectrum metadata `num`, `rt`, `tic` (when `metadata = true`)
followed by one column per m/z bin, named `mz_<centre>`. Being a `NamedTuple` of
vectors it is intrinsically a Tables source — `DataFrame(spectra_table(scans))`
works with no extra dependency on MassJ's side (mirrors [`isotope_table`](@ref)).

For a large number of bins, prefer the raw matrix from [`featurize`](@ref) (a
`NamedTuple` with thousands of fields is costly to build).

```julia
using DataFrames
df = DataFrame(spectra_table(scans; binsize = 1.0))
```
"""
function spectra_table(scans::AbstractVector{MSscans};
                       binsize::Real = 1.0,
                       mzrange::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
                       metadata::Bool = true)
    A = abundance_matrix(scans; binsize = binsize, threshold = 0.0,
                         mzrange = mzrange, drop_constant = false)
    X = A.matrix
    syms = Symbol[]
    cols = Any[]
    if metadata
        push!(syms, :num); push!(cols, [Int(_one(s.num, 0))     for s in scans])
        push!(syms, :rt);  push!(cols, [Float64(_one(s.rt, 0.0)) for s in scans])
        push!(syms, :tic); push!(cols, [Float64(s.tic)          for s in scans])
    end
    for (j, c) in enumerate(A.mz)
        push!(syms, Symbol("mz_", round(c, digits = 4)))
        push!(cols, X[:, j])
    end
    return NamedTuple{Tuple(syms)}(Tuple(cols))
end
