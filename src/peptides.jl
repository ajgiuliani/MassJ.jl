"""
Peptide fragment-ion calculator.

From an amino-acid sequence, compute fragment-ion m/z for the backbone series
a/b/c (N-terminal) and x/y/z (C-terminal), their ±H hydrogen-transfer ("radical")
variants common in UVPD / ExD, and the high-energy side-chain ions d, v, w.

Ion masses follow the Biemann / Roepstorff nomenclature as summarised by Mascot
(matrixscience.com): with [M] the residue-mass sum of the fragment and the
N-/C-terminal groups H / OH,

    a = [M] − CO       b = [M]            c = [M] + NH₃         (N-terminal, + H on the N-terminus)
    x = [M] + CO − H₂  y = [M] + H₂O      z = [M] + OH − NH₂    (C-terminal)

(written here in residue-sum form). Side-chain ions: `v` = y with the complete
side chain of the residue adjacent to the cleavage removed (that residue → Gly);
`d` = a and `w` = z with that residue's side chain cleaved at the Cβ–Cγ bond
(partial loss, modelled as homolytic so d/w are radicals). Cβ-branched residues
(Ile, Thr; Val degenerate) yield two isomers, d/d′ and w/w′.

Charge states are realised through the adduct machinery, so multiply-charged
fragments ([M+zH]ᶻ⁺) are electron-mass correct. Radical ETD/UVPD variants
(c·, z·, a+1 …) are requested through the `hshifts` keyword.
"""

export FragmentIon, fragment_ions, peptide_mass

# Residue elemental formulas (amino acid − H₂O) → monoisotopic masses on the
# shared element table, so fragment masses match `masses`/`adduct_mz` exactly.
const _AA_FORMULA = Dict{Char,String}(
    'G' => "C2H3NO",  'A' => "C3H5NO",   'S' => "C3H5NO2",  'P' => "C5H7NO",
    'V' => "C5H9NO",  'T' => "C4H7NO2",  'C' => "C3H5NOS",  'L' => "C6H11NO",
    'I' => "C6H11NO", 'N' => "C4H6N2O2", 'D' => "C4H5NO3",  'Q' => "C5H8N2O2",
    'K' => "C6H12N2O",'E' => "C5H7NO3",  'M' => "C5H9NOS",  'H' => "C6H7N3O",
    'F' => "C9H9NO",  'R' => "C6H12N4O", 'Y' => "C9H9NO2",  'W' => "C11H10N2O")

const _AA = Dict{Char,Float64}(c => masses(f)["Monoisotopic"] for (c, f) in _AA_FORMULA)

# Group masses on the same scale.
const _M_H   = first(Elements["H"]).m
const _M_CO  = masses("CO")["Monoisotopic"]
const _M_OH  = masses("OH")["Monoisotopic"]
const _M_NH2 = masses("NH2")["Monoisotopic"]
const _M_NH3 = masses("NH3")["Monoisotopic"]
const _M_H2O = masses("H2O")["Monoisotopic"]
const _M_CH3  = masses("CH3")["Monoisotopic"]
const _M_C2H5 = masses("C2H5")["Monoisotopic"]

# Complete side-chain loss (residue → Gly), used by v ions.
_v_loss(r::Char) = _AA[r] - _AA['G']

# Partial Cβ–Cγ side-chain losses, used by d and w ions. Returns the neutral
# loss mass(es): two for the Cβ-branched residues that give distinct isomers.
function _dw_losses(r::Char)
    r in ('G', 'A', 'P') && return Float64[]          # no Cβ–Cγ bond to cleave
    r == 'V' && return [_M_CH3]                        # two equivalent methyls
    r == 'I' && return [_M_CH3, _M_C2H5]               # d/d′
    r == 'T' && return [_M_CH3, _M_OH]                 # d/d′
    return [_AA[r] - (_AA['A'] - _M_H)]                # linear: side chain beyond CβH₂
end

"""
    struct FragmentIon

One computed peptide fragment ion.

    struct FragmentIon
        series::Symbol   # :a :b :c :x :y :z :d :v :w
        position::Int    # number of residues in the fragment
        charge::Int      # charge state z
        hshift::Int      # extra H atoms (radical / H-transfer); 0 = canonical
        isomer::Int      # 0 = primary, 1 = prime (d′/w′)
        mz::Float64
        label::String    # e.g. "b3", "y5(2+)", "c4+1", "d3", "w4'"
    end
"""
struct FragmentIon
    series::Symbol
    position::Int
    charge::Int
    hshift::Int
    isomer::Int
    mz::Float64
    label::String
end

function _label(series::Symbol, pos::Int, charge::Int, hshift::Int, isomer::Int)
    s = string(series) * string(pos)
    isomer == 1 && (s *= "'")
    hshift != 0 && (s *= hshift > 0 ? "+$hshift" : string(hshift))
    charge != 1 && (s *= "($(charge)+)")
    return s
end

"""
    peptide_mass(sequence) -> Float64

Neutral monoisotopic mass of a peptide (residue-mass sum + H₂O).
"""
function peptide_mass(sequence::AbstractString)
    seq = uppercase(strip(String(sequence)))
    m = _M_H2O
    for c in seq
        haskey(_AA, c) || error("peptide_mass: unknown residue '$c'.")
        m += _AA[c]
    end
    return m
end

# Emit one ion per (hshift, charge) for a given neutral fragment mass.
function _emit!(out, series, pos, neutral, charges, hshifts, adducts, isomer)
    for h in hshifts
        m = neutral + h * _M_H
        for z in charges
            mz = adduct_mz(m, adducts[z])
            push!(out, FragmentIon(series, pos, z, h, isomer, mz,
                                   _label(series, pos, z, h, isomer)))
        end
    end
end

const _NTERM = (:a, :b, :c, :d)
const _CTERM = (:x, :y, :z, :v, :w)

"""
    fragment_ions(sequence; ions=(:a,:b,:c,:x,:y,:z), charges=1:1, hshifts=(0,)) -> Vector{FragmentIon}

Compute fragment ions for `sequence` (one-letter amino-acid codes). `ions`
selects the series (add `:d`, `:v`, `:w` for the high-energy side-chain ions);
`charges` is the range of charge states; `hshifts` adds H-transfer variants
(e.g. `hshifts=(-1,0,1)` gives a−1/a/a+1, c−1/c/c+1, the z· radical, …).

# Examples
```julia-repl
julia> ions = fragment_ions("PEPTIDE"; ions=(:b,:y), charges=1:2);

julia> first(f for f in ions if f.label == "y3").mz
376.17144...
```
"""
function fragment_ions(sequence::AbstractString;
                       ions = (:a, :b, :c, :x, :y, :z),
                       charges = 1:1,
                       hshifts = (0,))
    seq = uppercase(strip(String(sequence)))
    res = collect(seq)
    n = length(res)
    n >= 2 || error("fragment_ions: sequence must have at least 2 residues.")
    for c in res
        haskey(_AA, c) || error("fragment_ions: unknown residue '$c' in sequence.")
    end
    rmass = [_AA[c] for c in res]

    adducts = Dict(z => Adduct("[M+$(z)H]$(z)+"; charge = z, add = "H$(z)") for z in charges)
    out = FragmentIon[]

    for ser in ions
        if ser in _NTERM
            for i in 1:n-1
                Nsum = sum(@view rmass[1:i])
                if ser === :d
                    for (k, L) in enumerate(_dw_losses(res[i]))
                        _emit!(out, :d, i, Nsum - _M_CO - L, charges, hshifts,
                               adducts, k == 2 ? 1 : 0)
                    end
                else
                    base = ser === :a ? Nsum - _M_CO :
                           ser === :b ? Nsum :
                                        Nsum + _M_NH3        # :c
                    _emit!(out, ser, i, base, charges, hshifts, adducts, 0)
                end
            end
        elseif ser in _CTERM
            for j in 1:n-1
                Csum = sum(@view rmass[n-j+1:n])
                r0   = res[n-j+1]                            # residue adjacent to cleavage
                if ser === :w
                    zbase = Csum + _M_OH - _M_NH2
                    for (k, L) in enumerate(_dw_losses(r0))
                        _emit!(out, :w, j, zbase - L, charges, hshifts, adducts,
                               k == 2 ? 1 : 0)
                    end
                elseif ser === :v
                    vloss = _v_loss(r0)
                    vloss == 0.0 && continue                 # Gly: v ≡ y, skip duplicate
                    _emit!(out, :v, j, Csum + _M_H2O - vloss, charges, hshifts, adducts, 0)
                else
                    base = ser === :x ? Csum + _M_OH + _M_CO - _M_H :
                           ser === :y ? Csum + _M_H2O :
                                        Csum + _M_OH - _M_NH2  # :z
                    _emit!(out, ser, j, base, charges, hshifts, adducts, 0)
                end
            end
        else
            error("fragment_ions: unknown ion series :$ser. Use :a :b :c :x :y :z :d :v :w.")
        end
    end
    return out
end
