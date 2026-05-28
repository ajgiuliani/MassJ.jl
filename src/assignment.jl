"""
Molecular-formula assignment from accurate mass, with isotope-pattern scoring.

The inverse of the forward isotope engine in isotopes.jl: enumerate candidate
formulas whose neutral monoisotopic mass matches an observed ion within a ppm
tolerance (bounded by element ranges and a ring-and-double-bond-equivalent
filter), then score each candidate's predicted isotope pattern against a measured
one. Reuses the adduct table (`neutral_mass`/`adduct_mz`) and `isotopic_distribution`.
"""

export FormulaCandidate, assign_formula, score_isotope_pattern

# Valences used for the ring + double-bond equivalents (RDBE) filter. Elements
# absent from this table contribute nothing (treated as valence 2).
const _VALENCE = Dict("C" => 4, "H" => 1, "N" => 3, "O" => 2, "S" => 2, "P" => 3,
                      "F" => 1, "Cl" => 1, "Br" => 1, "I" => 1, "Na" => 1, "K" => 1,
                      "Si" => 4, "B" => 3)

"""
    struct FormulaCandidate

A molecular formula proposed for an observed mass.

    struct FormulaCandidate
        formula::Dict{String,Int}   # element => count
        mass::Float64               # theoretical neutral monoisotopic mass
        error_ppm::Float64          # (theoretical − observed) / observed · 1e6
        rdbe::Float64               # ring + double-bond equivalents
    end
"""
struct FormulaCandidate
    formula::Dict{String,Int}
    mass::Float64
    error_ppm::Float64
    rdbe::Float64
end

# Ring + double-bond equivalents: 1 + ½·Σ nᵢ·(valenceᵢ − 2).
function _rdbe(f::Dict{String,Int})
    s = 0.0
    for (el, n) in f
        s += n * (get(_VALENCE, el, 2) - 2)
    end
    return 1.0 + 0.5 * s
end

# Branch-and-bound enumeration of element-count vectors whose monoisotopic mass
# lies within `tol_mass` of `target`. Kept as a top-level (non-closure) function
# so the recursion does not box captured state.
function _assign_recurse!(out::Vector{FormulaCandidate}, idx::Int, accmass::Float64,
                          counts::Vector{Int}, elsyms::Vector{String},
                          ranges::Vector{UnitRange{Int}}, emass::Vector{Float64},
                          suffixmax::Vector{Float64}, target::Float64, tol_mass::Float64,
                          rdbe::Tuple{Float64,Float64}, n::Int)
    if idx > n
        if abs(accmass - target) <= tol_mass
            f = Dict{String,Int}()
            @inbounds for k in 1:n
                counts[k] > 0 && (f[elsyms[k]] = counts[k])
            end
            d = _rdbe(f)
            if rdbe[1] <= d <= rdbe[2]
                push!(out, FormulaCandidate(f, accmass,
                                            (accmass - target) / target * 1e6, d))
            end
        end
        return
    end
    em = emass[idx]
    @inbounds for c in first(ranges[idx]):last(ranges[idx])
        m2 = accmass + c * em
        m2 - target > tol_mass && break                          # heavier from here on
        m2 + suffixmax[idx + 1] < target - tol_mass && continue   # too light even maxed
        counts[idx] = c
        _assign_recurse!(out, idx + 1, m2, counts, elsyms, ranges, emass,
                         suffixmax, target, tol_mass, rdbe, n)
    end
    counts[idx] = 0
    return
end

"""
    assign_formula(mz; adduct="[M+H]+", tol_ppm=5, elements=…, rdbe=(-0.5,40), max_results=50)

Propose molecular formulas for an ion observed at `mz`. The observed neutral
monoisotopic mass is recovered with [`neutral_mass`](@ref) for the given
`adduct`, then formulas within `tol_ppm` are enumerated over the supplied
`elements` (a `Dict` of element ⇒ count range, e.g. `"C" => 0:60`; a bare
integer `n` is read as `0:n`) and filtered by the ring + double-bond equivalents
range `rdbe`. Candidates are returned sorted by absolute mass error, truncated to
`max_results`.

# Examples
```julia-repl
julia> assign_formula(181.0707; adduct="[M+H]+", tol_ppm=5,
                      elements=Dict("C"=>0:12,"H"=>0:24,"O"=>0:12))[1].formula
Dict("C" => 6, "H" => 12, "O" => 6)
```
"""
function assign_formula(mz::Real;
                        adduct = "[M+H]+",
                        tol_ppm::Real = 5.0,
                        elements = Dict("C" => 0:60, "H" => 0:120, "N" => 0:10,
                                        "O" => 0:20, "P" => 0:5, "S" => 0:5),
                        rdbe::Tuple{<:Real,<:Real} = (-0.5, 40.0),
                        max_results::Int = 50)
    target = neutral_mass(mz, adduct)
    tol_mass = abs(target) * tol_ppm * 1e-6

    elsyms = String[]; ranges = UnitRange{Int}[]; emass = Float64[]
    for (el, r) in elements
        haskey(Elements, el) || error("Unknown element \"$el\" in `elements`.")
        push!(elsyms, el)
        push!(ranges, r isa Integer ? (0:Int(r)) : UnitRange{Int}(Int(first(r)), Int(last(r))))
        push!(emass, first(Elements[el]).m)
    end
    isempty(elsyms) && return FormulaCandidate[]

    # heaviest element first → tighter pruning
    o = sortperm(emass; rev = true)
    elsyms, ranges, emass = elsyms[o], ranges[o], emass[o]
    n = length(elsyms)

    suffixmax = zeros(Float64, n + 1)        # max mass reachable from elements > idx
    for i in n:-1:1
        suffixmax[i] = suffixmax[i + 1] + last(ranges[i]) * emass[i]
    end

    out = FormulaCandidate[]
    counts = zeros(Int, n)
    _assign_recurse!(out, 1, 0.0, counts, elsyms, ranges, emass, suffixmax,
                     target, tol_mass, (Float64(rdbe[1]), Float64(rdbe[2])), n)

    sort!(out, by = c -> abs(c.error_ppm))
    return length(out) > max_results ? out[1:max_results] : out
end

"""
    score_isotope_pattern(mz, intensity, candidate; adduct="[M+H]+", tol=0.02, p_target=0.99)
    score_isotope_pattern(spec::MSscans, candidate; …)

Cosine similarity in `[0, 1]` between a measured isotope pattern (`mz` /
`intensity`, or the peaks of an [`MSscans`](@ref)) and the pattern predicted for
`candidate`. The candidate's isotopologues are produced with
[`isotopic_distribution`](@ref) and placed at observed m/z through the `adduct`;
each predicted peak is matched to the nearest measured peak within `tol` (Da).
"""
function score_isotope_pattern(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real},
                               candidate::FormulaCandidate;
                               adduct = "[M+H]+", tol::Real = 0.02, p_target::Real = 0.99)
    (isempty(mz) || isempty(candidate.formula)) && return 0.0
    theo = isotopic_distribution(candidate.formula, p_target; charge = 1)
    # row 1 holds labels; data rows: col 1 = neutral isotopologue mass, col 2 = probability
    masses_iso = Float64.(@view theo[2:end, 1])
    probs      = Float64.(@view theo[2:end, 2])

    vt = Float64[]; vm = Float64[]
    for (m, p) in zip(masses_iso, probs)
        tmz = adduct_mz(m, adduct)
        j = num2pnt(mz, tmz)
        push!(vt, p)
        push!(vm, abs(mz[j] - tmz) <= tol ? float(intensity[j]) : 0.0)
    end
    (isempty(vt) || all(iszero, vm)) && return 0.0
    return dot(vt, vm) / (norm(vt) * norm(vm))
end

score_isotope_pattern(spec::MSscans, candidate::FormulaCandidate; kwargs...) =
    score_isotope_pattern(spec.mz, spec.int, candidate; kwargs...)
