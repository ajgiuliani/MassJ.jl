"""
Single-pass filtering via composed predicates.

Each `FilterType` is converted to a predicate function `scan::MSscans -> Bool` via
`to_predicate`. Multiple predicates are composed into a single function that
short-circuits on the first `false`, enabling single-pass filter + accumulate.

RT-based filters that need access to the full retention time array are handled by
`to_predicate(scans, filter)` which pre-computes index bounds once.
"""


# --- Generic predicate composition ---

"""
    to_predicate(f::FilterType)
Convert a `FilterType` to a predicate `scan::MSscans -> Bool`.
Falls back to the two-argument form when the filter needs context from the full scan list.
"""
function to_predicate end

"""
    to_predicate(scans::AbstractVector{MSscans}, f::FilterType)
Convert a `FilterType` to a predicate, with access to the full scan list for
filters that need global context (e.g. RT needs the retention time array).
Defaults to ignoring the scans argument.
"""
to_predicate(scans::AbstractVector{MSscans}, f::FilterType) = to_predicate(f)

"""
    compose_predicates(scans::AbstractVector{MSscans}, filters::Tuple{Vararg{FilterType}})
Build a single predicate from multiple `FilterType`s. Returns `scan -> Bool`.
"""
function compose_predicates(scans::AbstractVector{MSscans}, filters)
    isempty(filters) && return _ -> true
    preds = Tuple(to_predicate(scans, f) for f in filters)
    return scan -> all(p -> p(scan), preds)
end


# Each predicate inspects a single scan — an `MSscans` with length-1 provenance —
# so the per-scan value is the first (and only) element of each provenance vector.

# --- Level ---

to_predicate(f::Level{<:Int}) = scan -> scan.level[1] == f.arg
to_predicate(f::Level{<:AbstractVector}) = scan -> scan.level[1] ∈ f.arg


# --- Precursor ---

to_predicate(f::Precursor{<:Real}) = scan -> scan.precursor[1] == f.arg
to_predicate(f::Precursor{<:AbstractVector}) = scan -> scan.precursor[1] ∈ f.arg


# --- Activation_Energy ---

to_predicate(f::Activation_Energy{<:Real}) = scan -> scan.collisionEnergy[1] == f.arg
to_predicate(f::Activation_Energy{<:AbstractVector}) = scan -> scan.collisionEnergy[1] ∈ f.arg


# --- Activation_Method ---

to_predicate(f::Activation_Method{<:String}) = scan -> scan.activationMethod[1] == f.arg
to_predicate(f::Activation_Method{<:AbstractVector}) = scan -> scan.activationMethod[1] ∈ f.arg


# --- Polarity ---

to_predicate(f::Polarity{<:String}) = scan -> scan.polarity[1] == f.arg
to_predicate(f::Polarity{<:AbstractVector}) = scan -> scan.polarity[1] ∈ f.arg


# --- Scan ---

to_predicate(f::Scan{<:Int}) = scan -> scan.num[1] == f.arg
to_predicate(f::Scan{<:AbstractVector}) = scan -> scan.num[1] ∈ f.arg


# --- IC (ion current range) ---

to_predicate(f::IC{<:AbstractVector}) = scan -> f.arg[1] <= scan.tic <= f.arg[2]


# --- DriftTime ---

to_predicate(f::DriftTime{<:Real}) = scan -> scan.driftTime[1] == f.arg
to_predicate(f::DriftTime{<:AbstractVector}) = scan -> f.arg[1] <= scan.driftTime[1] <= f.arg[2]


# --- CompensationVoltage ---

to_predicate(f::CompensationVoltage{<:Real}) = scan -> scan.compensationVoltage[1] == f.arg
to_predicate(f::CompensationVoltage{<:AbstractVector}) = scan -> f.arg[1] <= scan.compensationVoltage[1] <= f.arg[2]


# --- ChargeState ---

to_predicate(f::ChargeState{<:Int}) = scan -> scan.chargeState[1] == f.arg
to_predicate(f::ChargeState{<:AbstractVector}) = scan -> scan.chargeState[1] ∈ f.arg


# --- SpectrumType (scalar field) ---

to_predicate(f::SpectrumType{Symbol}) = scan -> scan.spectrumType == f.arg
to_predicate(f::SpectrumType{<:AbstractVector}) = scan -> scan.spectrumType ∈ f.arg


# --- MobilityType (scalar field) ---

to_predicate(f::MobilityType{Symbol}) = scan -> scan.mobilityType == f.arg
to_predicate(f::MobilityType{<:AbstractVector}) = scan -> scan.mobilityType ∈ f.arg


# --- MetaData (generic per-spectrum cvParam) ---

# presence test
to_predicate(f::MetaData{Nothing}) = scan -> haskey(scan.metadata, f.key)
# numeric exact match
to_predicate(f::MetaData{<:Real}) =
    scan -> haskey(scan.metadata, f.key) &&
            (v = scan.metadata[f.key]; v isa Real && v == f.value)
# numeric range [lo, hi]
to_predicate(f::MetaData{<:AbstractVector{<:Real}}) =
    scan -> haskey(scan.metadata, f.key) &&
            (v = scan.metadata[f.key]; v isa Real && f.value[1] <= v <= f.value[2])
# substring match on a stored string
to_predicate(f::MetaData{<:AbstractString}) =
    scan -> haskey(scan.metadata, f.key) && occursin(f.value, string(scan.metadata[f.key]))


# --- RT (needs global context: retention time array → index mapping) ---

function to_predicate(scans::AbstractVector{MSscans}, f::RT{<:Real})
    rt = retention_time(scans)
    target_num = num2pnt(rt, f.arg)
    return scan -> scan.num[1] == target_num
end

function to_predicate(scans::AbstractVector{MSscans}, f::RT{<:AbstractVector{<:Real}})
    rt = retention_time(scans)
    bounds = Set{Int}()
    for i in 1:2:length(f.arg)
        lo = num2pnt(rt, f.arg[i])
        hi = num2pnt(rt, f.arg[i+1])
        for idx in lo:hi
            push!(bounds, idx)
        end
    end
    return scan -> scan.num[1] ∈ bounds
end

function to_predicate(scans::AbstractVector{MSscans}, f::RT{<:AbstractVector{<:AbstractVector}})
    rt = retention_time(scans)
    bounds = Set{Int}()
    for el in f.arg
        lo = num2pnt(rt, el[1])
        hi = num2pnt(rt, el[2])
        for idx in lo:hi
            push!(bounds, idx)
        end
    end
    return scan -> scan.num[1] ∈ bounds
end


# --- InstrumentConfig (needs global context: the run-level instrument table) ---

# True if any cvParam dict in `cvs` has a name/accession/value equal to `query`.
_cvparams_match(cvs, query::AbstractString) =
    any(cv -> get(cv, "name", "") == query ||
              get(cv, "accession", "") == query ||
              string(get(cv, "value", "")) == query, cvs)

# Set of instrument-configuration ids matching `query` (by id, cvParam, or component
# type), resolved from an MSrun's instrument table; `nothing` when no table is available.
function _instrument_config_matches(scans, query::AbstractString)
    (scans isa MSrun && haskey(scans.metadata, "instruments")) || return nothing
    ids = Set{String}()
    for ic in scans.metadata["instruments"]
        id = get(ic, "id", "")
        matched = id == query || _cvparams_match(get(ic, "cv_params", Any[]), query)
        if !matched
            for comp in get(ic, "components", Any[])
                if get(comp, "type", "") == query || _cvparams_match(get(comp, "cv_params", Any[]), query)
                    matched = true
                    break
                end
            end
        end
        matched && push!(ids, id)
    end
    return ids
end

function to_predicate(scans::AbstractVector{MSscans}, f::InstrumentConfig)
    matches = _instrument_config_matches(scans, f.query)
    if matches === nothing
        # No run-level table: match the per-spectrum ref directly against the query.
        return scan -> get(scan.metadata, "instrument_configuration_ref", "") == f.query
    end
    return scan -> get(scan.metadata, "instrument_configuration_ref", "") ∈ matches
end
