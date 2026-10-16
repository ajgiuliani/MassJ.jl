# ============================================================================
# Universal Spectrum Identifier (USI) retrieval
#
# Resolve a single spectrum referenced by a USI (PSI standard,
# https://www.psidev.info/usi) from a PROXI-compliant server and return it as
# an `MSscans`. This lets a spectrum be loaded by a portable reference --- the
# "permalink" of a spectrum in a public repository --- rather than from a local
# file, mirroring `spectrum_utils.MsmsSpectrum.from_usi`.
# ============================================================================

import Downloads
import JSON

export load_usi

"""Default PROXI aggregator that dispatches a USI to the repository hosting it."""
const USI_DEFAULT_RESOLVER =
    "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra"

# percent-encode a (ASCII) USI for safe inclusion in a query string
function _percent_encode(s::AbstractString)
    io = IOBuffer()
    for b in codeunits(s)
        c = Char(b)
        if isascii(c) && (isletter(c) || isdigit(c) || c in "-._~")
            write(io, c)
        else
            write(io, '%', uppercase(string(b, base = 16, pad = 2)))
        end
    end
    String(take!(io))
end

# value coercion (PROXI attribute values may be strings or numbers)
_to_f(::Nothing) = nothing
_to_f(x::Real) = Float64(x)
_to_f(x::AbstractString) = tryparse(Float64, x)
_to_i(::Nothing) = nothing
_to_i(x::Real) = round(Int, x)
_to_i(x::AbstractString) = tryparse(Int, x)

# look up a PROXI attribute value by cvParam accession (or `nothing`)
function _usi_attr(attrs, accession)
    for a in attrs
        a isa AbstractDict && get(a, "accession", "") == accession &&
            return get(a, "value", nothing)
    end
    return nothing
end

# Thermo-style filter string fallback, e.g.
# "FTMS + p NSI d Full ms2 767.97@hcd32.00 [110.00-1550.00]"
function _parse_filter_string(fs::AbstractString)
    level = activation = ce = polarity = nothing
    m = match(r"ms(\d+)"i, fs);                 m === nothing || (level = parse(Int, m.captures[1]))
    m = match(r"@([a-z]+)(\d+(?:\.\d+)?)"i, fs)
    if m !== nothing
        activation = uppercase(m.captures[1]); ce = parse(Float64, m.captures[2])
    end
    occursin(r"\s\+\s", fs) && (polarity = "+")
    occursin(r"\s-\s", fs)  && (polarity = "-")
    return (level = level, activation = activation, ce = ce, polarity = polarity)
end

# explicit activation cvParam (PROXI sometimes carries one)
const _USI_ACTIVATION = Dict("MS:1000133" => "CID", "MS:1000422" => "HCD",
                             "MS:1000598" => "ETD", "MS:1000250" => "ECD",
                             "MS:1000599" => "PQD", "MS:1000262" => "IRMPD",
                             "MS:1003246" => "UVPD")
function _activation_from_attrs(attrs)
    for a in attrs
        a isa AbstractDict || continue
        acc = get(a, "accession", "")
        haskey(_USI_ACTIVATION, acc) && return _USI_ACTIVATION[acc]
    end
    return nothing
end

function _http_get(url::AbstractString)
    buf = IOBuffer()
    try
        Downloads.download(url, buf; headers = ["Accept" => "application/json"])
    catch e
        error("USI retrieval failed (network or resolver error): " *
              sprint(showerror, e))
    end
    return String(take!(buf))
end

# Build an `MSscans` from one PROXI spectrum object (Dict with mzs/intensities/attributes).
function _msscans_from_proxi(spec::AbstractDict, usi::AbstractString)
    mzs   = Float64.(spec["mzs"])
    ints  = Float64.(spec["intensities"])
    attrs = get(spec, "attributes", Any[])
    A(acc) = _usi_attr(attrs, acc)

    fs = A("MS:1000512")
    fl = fs === nothing ? (level = nothing, activation = nothing, ce = nothing, polarity = nothing) :
                          _parse_filter_string(String(fs))

    prec    = something(_to_f(A("MS:1000744")), 0.0)                 # selected ion m/z
    charge  = something(_to_i(A("MS:1000041")), 0)                  # charge state
    level   = something(_to_i(A("MS:1000511")), fl.level, prec > 0 ? 2 : 1)
    polar   = something(A("MS:1000130") === nothing ? nothing : "+",
                        A("MS:1000129") === nothing ? nothing : "-",
                        fl.polarity, "+")
    activ   = something(_activation_from_attrs(attrs), fl.activation, "")
    ce      = something(_to_f(A("MS:1000045")), fl.ce, 0.0)         # collision energy
    rt_s    = _to_f(A("MS:1000016"))                               # scan start time (s)
    rt_min  = rt_s === nothing ? 0.0 : rt_s / 60.0                 # MassJ stores rt in minutes
    scanno  = something(_to_i(A("MS:1008025")), _to_i(A("MS:1003057")), 1)

    tic = isempty(ints) ? 0.0 : sum(ints)
    bi  = isempty(ints) ? 0 : argmax(ints)
    bpmz, bpint = bi == 0 ? (0.0, 0.0) : (mzs[bi], ints[bi])

    md = Dict{String,Any}("usi" => String(usi), "source_format" => "USI")
    fs === nothing || (md["filter_string"] = String(fs))
    for (key, acc) in (("isolation_window_target_mz", "MS:1000827"),
                       ("isolation_window_lower_offset", "MS:1000828"),
                       ("isolation_window_upper_offset", "MS:1000829"))
        v = _to_f(A(acc)); v === nothing || (md[key] = v)
    end

    return MSscans(scanno, rt_min, tic, mzs, ints, level,
                   bpmz, bpint, prec, polar, activ, ce,
                   charge, :centroid, -1.0, 0.0, :none, md)
end

"""
    load_usi(usi; resolver = USI_DEFAULT_RESOLVER, verbose = false) -> MSscans

Retrieve the single mass spectrum referenced by a Universal Spectrum Identifier
(USI) from a PROXI-compliant server and return it as an [`MSscans`](@ref).

A USI is a portable reference to one spectrum in a public proteomics or
metabolomics repository, e.g.

    mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555

The peak list, precursor m/z, charge, activation, collision energy, and isolation
window (when provided by the server) are populated; the originating USI and any
instrument filter string are kept in `metadata`. The default `resolver` is the
ProteomeCentral aggregator, which dispatches the USI to the hosting repository;
pass another PROXI endpoint (e.g. a MassIVE or PeptideAtlas URL) to override it.

```julia
spec = load_usi("mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555")
spec.precursor[1]    # precursor m/z
spec.mz, spec.int    # fragment peaks
```

Network access is required. Throws if the USI is malformed, the resolver is
unreachable, or the spectrum is not found.
"""
function load_usi(usi::AbstractString;
                  resolver::AbstractString = USI_DEFAULT_RESOLVER,
                  verbose::Bool = false)
    startswith(usi, "mzspec:") ||
        error("not a USI (must start with \"mzspec:\"): $usi")
    url = string(resolver, "?usi=", _percent_encode(String(usi)), "&resultType=full")
    verbose && @info "Resolving USI" url
    data = JSON.parse(_http_get(url))

    if data isa AbstractDict                       # error object, not a spectrum list
        msg = something(get(data, "detail", nothing), get(data, "title", nothing),
                        get(data, "message", nothing), "unknown error")
        error("PROXI resolver returned an error for USI $usi: $msg")
    end
    data isa AbstractVector ||
        error("unexpected PROXI response for USI $usi")

    for s in data
        s isa AbstractDict && !isempty(get(s, "mzs", Any[])) &&
            return _msscans_from_proxi(s, usi)
    end
    error("USI not found, or resolver returned no peaks: $usi")
end
