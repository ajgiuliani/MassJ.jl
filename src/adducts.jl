"""
Adduct and neutral-loss handling.

Converts between neutral monoisotopic masses and observed ion m/z for a set of
common electrospray adducts, reusing `formula` and `masses` from isotopes.jl to
turn the added/removed groups into mass shifts.
"""

export Adduct, Adducts, adduct_mz, neutral_mass

# Electron rest mass (Da), CODATA. Accounting for it matters at the sub-ppm
# level; the legacy `deconv` adduct convention ignores it (see `adduct_mz`).
const m_electron = 0.000548579909065

"""
    struct Adduct

Description of a mass-spectrometric adduct ion. The ion mass is

    ion_mass = n·M + mass(add) − mass(sub) − charge·mₑ

where `M` is the analyte neutral monoisotopic mass and `mₑ` the electron mass,
and the observed value is `m/z = ion_mass / |charge|`.

    struct Adduct
        name::String     # display name, e.g. "[M+H]+"
        charge::Int      # signed charge (e.g. +1, +2, -1)
        n::Int           # analyte multiplicity (1 = monomer, 2 = dimer, ...)
        add::String      # formula added ("" if none)
        sub::String      # formula removed ("" if none)
    end
"""
struct Adduct
    name::String
    charge::Int
    n::Int
    add::String
    sub::String
    Adduct(name::AbstractString, charge::Integer, n::Integer,
           add::AbstractString, sub::AbstractString) =
        new(String(name), Int(charge), Int(n), String(add), String(sub))
end

"""
    Adduct(name; charge, n=1, add="", sub="")

Keyword constructor for an [`Adduct`](@ref).
"""
Adduct(name::AbstractString; charge::Integer, n::Integer = 1,
       add::AbstractString = "", sub::AbstractString = "") =
    Adduct(name, charge, n, add, sub)

# Monoisotopic mass of a (possibly empty) formula fragment.
_fragment_mass(f::AbstractString) = isempty(f) ? 0.0 : masses(formula(f))["Monoisotopic"]

"""
    polarity(a::Adduct)

Return `"+"` or `"-"` for the sign of the adduct charge.
"""
polarity(a::Adduct) = a.charge < 0 ? "-" : "+"

"""
    Adducts

Lookup table of common adducts keyed by name. Use [`adduct_mz`](@ref) and
[`neutral_mass`](@ref) with a key from this table or a custom [`Adduct`](@ref).
"""
const Adducts = Dict{String,Adduct}(
    # positive mode
    "[M+H]+"      => Adduct("[M+H]+",      charge = +1, add = "H"),
    "[M+2H]2+"    => Adduct("[M+2H]2+",    charge = +2, add = "H2"),
    "[M+3H]3+"    => Adduct("[M+3H]3+",    charge = +3, add = "H3"),
    "[M+Na]+"     => Adduct("[M+Na]+",     charge = +1, add = "Na"),
    "[M+K]+"      => Adduct("[M+K]+",      charge = +1, add = "K"),
    "[M+NH4]+"    => Adduct("[M+NH4]+",    charge = +1, add = "NH4"),
    "[M+H-H2O]+"  => Adduct("[M+H-H2O]+",  charge = +1, add = "H",   sub = "H2O"),
    "[M+2Na-H]+"  => Adduct("[M+2Na-H]+",  charge = +1, add = "Na2", sub = "H"),
    "[2M+H]+"     => Adduct("[2M+H]+",     charge = +1, n = 2, add = "H"),
    "[2M+Na]+"    => Adduct("[2M+Na]+",    charge = +1, n = 2, add = "Na"),
    "[M]+"        => Adduct("[M]+",        charge = +1),
    # negative mode
    "[M-H]-"      => Adduct("[M-H]-",      charge = -1, sub = "H"),
    "[M-2H]2-"    => Adduct("[M-2H]2-",    charge = -2, sub = "H2"),
    "[M-3H]3-"    => Adduct("[M-3H]3-",    charge = -3, sub = "H3"),
    "[M+Cl]-"     => Adduct("[M+Cl]-",     charge = -1, add = "Cl"),
    "[M+CHO2]-"   => Adduct("[M+CHO2]-",   charge = -1, add = "CH2O2", sub = "H"),  # formate
    "[M+C2H3O2]-" => Adduct("[M+C2H3O2]-", charge = -1, add = "C2H4O2", sub = "H"), # acetate
    "[M-H-H2O]-"  => Adduct("[M-H-H2O]-",  charge = -1, sub = "H3O"),
    "[2M-H]-"     => Adduct("[2M-H]-",     charge = -1, n = 2, sub = "H"),
    "[M]-"        => Adduct("[M]-",        charge = -1),
)

# Resolve a name (looked up in `Adducts`) or pass through an `Adduct`.
_resolve_adduct(a::Adduct) = a
function _resolve_adduct(name::AbstractString)
    haskey(Adducts, name) && return Adducts[name]
    error("Unknown adduct \"$name\". Known adducts: $(join(sort(collect(keys(Adducts))), ", ")). " *
          "Build a custom one with Adduct(name; charge, n, add, sub).")
end

"""
    adduct_mz(neutral_mass, adduct; electron=true)

Return the observed `m/z` for an analyte of monoisotopic mass `neutral_mass`
ionised as `adduct` (an [`Adduct`](@ref) or a name from [`Adducts`](@ref)). When
`electron=false`, the electron mass is ignored (matching the legacy `deconv`
convention).

# Examples
```julia-repl
julia> adduct_mz(180.06339, "[M+H]+")        # glucose, protonated
181.07067...

julia> adduct_mz(180.06339, "[M+Na]+")
203.05261...
```
"""
function adduct_mz(neutral_mass::Real, adduct; electron::Bool = true)
    a = _resolve_adduct(adduct)
    me = electron ? m_electron : 0.0
    ion = a.n * neutral_mass + _fragment_mass(a.add) - _fragment_mass(a.sub) - a.charge * me
    return ion / abs(a.charge)
end

"""
    neutral_mass(mz, adduct; electron=true)

Inverse of [`adduct_mz`](@ref): return the analyte neutral monoisotopic mass
that would be observed at `mz` for the given `adduct`.

# Examples
```julia-repl
julia> neutral_mass(181.07067, "[M+H]+")
180.06339...
```
"""
function neutral_mass(mz::Real, adduct; electron::Bool = true)
    a = _resolve_adduct(adduct)
    me = electron ? m_electron : 0.0
    ion = mz * abs(a.charge)
    return (ion - _fragment_mass(a.add) + _fragment_mass(a.sub) + a.charge * me) / a.n
end
