"""
Single-pass filtering via composed predicates.

Each `FilterType` is converted to a predicate function `scan::MSscan -> Bool` via
`to_predicate`. Multiple predicates are composed into a single function that
short-circuits on the first `false`, enabling single-pass filter + accumulate.

RT-based filters that need access to the full retention time array are handled by
`to_predicate(scans, filter)` which pre-computes index bounds once.
"""


# --- Generic predicate composition ---

"""
    to_predicate(f::FilterType)
Convert a `FilterType` to a predicate `scan::MSscan -> Bool`.
Falls back to the two-argument form when the filter needs context from the full scan list.
"""
function to_predicate end

"""
    to_predicate(scans::Vector{MSscan}, f::FilterType)
Convert a `FilterType` to a predicate, with access to the full scan list for
filters that need global context (e.g. RT needs the retention time array).
Defaults to ignoring the scans argument.
"""
to_predicate(scans::Vector{MSscan}, f::FilterType) = to_predicate(f)

"""
    compose_predicates(scans::Vector{MSscan}, filters::Tuple{Vararg{FilterType}})
Build a single predicate from multiple `FilterType`s. Returns `scan -> Bool`.
"""
function compose_predicates(scans::Vector{MSscan}, filters)
    isempty(filters) && return _ -> true
    preds = Tuple(to_predicate(scans, f) for f in filters)
    return scan -> all(p -> p(scan), preds)
end


# --- Level ---

to_predicate(f::Level{<:Int}) = scan -> scan.level == f.arg
to_predicate(f::Level{<:AbstractVector}) = scan -> scan.level ∈ f.arg


# --- Precursor ---

to_predicate(f::Precursor{<:Real}) = scan -> scan.precursor == f.arg
to_predicate(f::Precursor{<:AbstractVector}) = scan -> scan.precursor ∈ f.arg


# --- Activation_Energy ---

to_predicate(f::Activation_Energy{<:Real}) = scan -> scan.collisionEnergy == f.arg
to_predicate(f::Activation_Energy{<:AbstractVector}) = scan -> scan.collisionEnergy ∈ f.arg


# --- Activation_Method ---

to_predicate(f::Activation_Method{<:String}) = scan -> scan.activationMethod == f.arg
to_predicate(f::Activation_Method{<:AbstractVector}) = scan -> scan.activationMethod ∈ f.arg


# --- Polarity ---

to_predicate(f::Polarity{<:String}) = scan -> scan.polarity == f.arg
to_predicate(f::Polarity{<:AbstractVector}) = scan -> scan.polarity ∈ f.arg


# --- Scan ---

to_predicate(f::Scan{<:Int}) = scan -> scan.num == f.arg
to_predicate(f::Scan{<:AbstractVector}) = scan -> scan.num ∈ f.arg


# --- IC (ion current range) ---

to_predicate(f::IC{<:AbstractVector}) = scan -> f.arg[1] <= scan.tic <= f.arg[2]


# --- DriftTime ---

to_predicate(f::DriftTime{<:Real}) = scan -> scan.driftTime == f.arg
to_predicate(f::DriftTime{<:AbstractVector}) = scan -> f.arg[1] <= scan.driftTime <= f.arg[2]


# --- CompensationVoltage ---

to_predicate(f::CompensationVoltage{<:Real}) = scan -> scan.compensationVoltage == f.arg
to_predicate(f::CompensationVoltage{<:AbstractVector}) = scan -> f.arg[1] <= scan.compensationVoltage <= f.arg[2]


# --- RT (needs global context: retention time array → index mapping) ---

function to_predicate(scans::Vector{MSscan}, f::RT{<:Real})
    rt = retention_time(scans)
    target_num = num2pnt(rt, f.arg)
    return scan -> scan.num == target_num
end

function to_predicate(scans::Vector{MSscan}, f::RT{<:AbstractVector{<:Real}})
    rt = retention_time(scans)
    bounds = Set{Int}()
    for i in 1:2:length(f.arg)
        lo = num2pnt(rt, f.arg[i])
        hi = num2pnt(rt, f.arg[i+1])
        for idx in lo:hi
            push!(bounds, idx)
        end
    end
    return scan -> scan.num ∈ bounds
end

function to_predicate(scans::Vector{MSscan}, f::RT{<:AbstractVector{<:AbstractVector}})
    rt = retention_time(scans)
    bounds = Set{Int}()
    for el in f.arg
        lo = num2pnt(rt, el[1])
        hi = num2pnt(rt, el[2])
        for idx in lo:hi
            push!(bounds, idx)
        end
    end
    return scan -> scan.num ∈ bounds
end
