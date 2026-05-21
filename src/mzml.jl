"""
Interface to the mzML file format (PSI standard).
Uses LightXML for parsing and the same binary decoding pipeline as mzXML (Codecs, Libz).
mzML stores binary data in little-endian byte order (unlike mzXML which uses network/big-endian).
"""


# CV accession constants for lookup
const CV_MS_LEVEL          = "MS:1000511"
const CV_POSITIVE_SCAN     = "MS:1000130"
const CV_NEGATIVE_SCAN     = "MS:1000129"
const CV_CENTROID          = "MS:1000127"
const CV_PROFILE           = "MS:1000128"
const CV_TIC               = "MS:1000285"
const CV_BASE_PEAK_MZ      = "MS:1000504"
const CV_BASE_PEAK_INT     = "MS:1000505"
const CV_SCAN_START_TIME   = "MS:1000016"
const CV_SELECTED_ION_MZ   = "MS:1000744"
const CV_CHARGE_STATE      = "MS:1000041"
const CV_COLLISION_ENERGY  = "MS:1000045"
const CV_MZ_ARRAY          = "MS:1000514"
const CV_INT_ARRAY         = "MS:1000515"
const CV_64BIT             = "MS:1000523"
const CV_32BIT             = "MS:1000521"
const CV_ZLIB              = "MS:1000574"
const CV_NO_COMPRESSION    = "MS:1000576"
const CV_UNIT_MINUTE       = "UO:0000031"
const CV_UNIT_SECOND       = "UO:0000010"
const CV_DRIFT_TIME        = "MS:1002476"
const CV_INV_K0            = "MS:1002815"
const CV_FAIMS_CV          = "MS:1001581"
const CV_SELEXION_CV       = "MS:1003371"

# Activation method accessions
const CV_CID               = "MS:1000133"
const CV_HCD               = "MS:1000422"
const CV_ETD               = "MS:1000598"
const CV_ECD               = "MS:1000250"
const CV_PQD               = "MS:1000599"
const CV_IRMPD             = "MS:1000262"
const CV_SID               = "MS:1000282"
const CV_UVPD              = "MS:1003246"

const ACTIVATION_METHODS = Dict(
    CV_CID   => "CID",
    CV_HCD   => "HCD",
    CV_ETD   => "ETD",
    CV_ECD   => "ECD",
    CV_PQD   => "PQD",
    CV_IRMPD => "IRMPD",
    CV_SID   => "SID",
    CV_UVPD  => "UVPD",
)


"""
    get_cv_param(elem::XMLElement, accession::String)
Search for a cvParam child with the given accession. Returns the XMLElement or nothing.
"""
function get_cv_param(elem::XMLElement, accession::String)
    for child in child_elements(elem)
        if name(child) == "cvParam" && attribute(child, "accession") == accession
            return child
        end
    end
    return nothing
end

"""
    get_cv_value(elem::XMLElement, accession::String, default="")
Get the value attribute of a cvParam with the given accession. Returns default if not found.
"""
function get_cv_value(elem::XMLElement, accession::String, default::String="")
    cv = get_cv_param(elem, accession)
    if cv !== nothing
        val = attribute(cv, "value")
        return val !== nothing ? val : default
    end
    return default
end

"""
    has_cv_param(elem::XMLElement, accession::String)
Check whether elem has a cvParam child with the given accession.
"""
function has_cv_param(elem::XMLElement, accession::String)
    return get_cv_param(elem, accession) !== nothing
end


"""
    find_mzml_root(xdoc::XMLDocument)
Navigate to the mzML element, handling both raw mzML and indexedmzML wrappers.
Returns the mzML element.
"""
function find_mzml_root(xdoc::XMLDocument)
    xroot = root(xdoc)
    rname = name(xroot)
    if rname == "mzML"
        return xroot
    elseif rname == "indexedmzML"
        mzml = find_element(xroot, "mzML")
        if mzml !== nothing
            return mzml
        end
    end
    error("Not an mzML file.")
end


"""
    info_mzml(filename::String, info::Vector{String}, verbose::Bool=false)
Returns the information content of an mzML file into a string.
"""
function info_mzml(filename::String, info::Vector{String}, verbose::Bool=false)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)

    if verbose
        # File description
        fileDesc = find_element(mzml, "fileDescription")
        if fileDesc !== nothing
            srcList = find_element(fileDesc, "sourceFileList")
            if srcList !== nothing
                for sf in child_elements(srcList)
                    if name(sf) == "sourceFile"
                        fname = attribute(sf, "name")
                        if fname !== nothing
                            push!(info, "parentFile: " * fname)
                        end
                    end
                end
            end
        end

        # Instrument configuration
        instrList = find_element(mzml, "instrumentConfigurationList")
        if instrList !== nothing
            for ic in child_elements(instrList)
                if name(ic) == "instrumentConfiguration"
                    for cv in child_elements(ic)
                        if name(cv) == "cvParam"
                            n = attribute(cv, "name")
                            if n !== nothing
                                push!(info, "msModel: " * n)
                            end
                        end
                    end
                    compList = find_element(ic, "componentList")
                    if compList !== nothing
                        for comp in child_elements(compList)
                            cname = name(comp)
                            for cv in child_elements(comp)
                                if name(cv) == "cvParam"
                                    n = attribute(cv, "name")
                                    if n !== nothing
                                        if cname == "source"
                                            push!(info, "msIonisation: " * n)
                                        elseif cname == "analyzer"
                                            push!(info, "msMassAnalyzer: " * n)
                                        elseif cname == "detector"
                                            push!(info, "msDetector: " * n)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        # Software
        swList = find_element(mzml, "softwareList")
        if swList !== nothing
            for sw in child_elements(swList)
                if name(sw) == "software"
                    sid = attribute(sw, "id")
                    ver = attribute(sw, "version")
                    sid = sid !== nothing ? sid : "unknown"
                    ver = ver !== nothing ? ver : ""
                    push!(info, "software: " * sid * " " * ver)
                end
            end
        end
    end

    # Count spectra and collect unique scan descriptions
    run_elem = find_element(mzml, "run")
    if run_elem === nothing
        free(xdoc)
        return info
    end

    specList = find_element(run_elem, "spectrumList")
    if specList === nothing
        free(xdoc)
        return info
    end

    scanCount = attribute(specList, "count")
    if scanCount !== nothing
        push!(info, scanCount * " scans")
    end

    seen = Set{String}()
    for spec in child_elements(specList)
        if name(spec) != "spectrum"
            continue
        end

        msLevel = get_cv_value(spec, CV_MS_LEVEL, "0")
        filter = "MS" * msLevel

        if has_cv_param(spec, CV_POSITIVE_SCAN)
            filter *= "+"
        elseif has_cv_param(spec, CV_NEGATIVE_SCAN)
            filter *= "-"
        end

        # Precursor info
        precList = find_element(spec, "precursorList")
        if precList !== nothing
            for prec in child_elements(precList)
                if name(prec) != "precursor"
                    continue
                end
                selIonList = find_element(prec, "selectedIonList")
                if selIonList !== nothing
                    for selIon in child_elements(selIonList)
                        if name(selIon) == "selectedIon"
                            pmz = get_cv_value(selIon, CV_SELECTED_ION_MZ)
                            if pmz != ""
                                filter *= " " * pmz
                            end
                        end
                    end
                end
                act = find_element(prec, "activation")
                if act !== nothing
                    actMethod = ""
                    ce = ""
                    for cv in child_elements(act)
                        if name(cv) == "cvParam"
                            acc = attribute(cv, "accession")
                            if acc !== nothing && haskey(ACTIVATION_METHODS, acc)
                                actMethod = ACTIVATION_METHODS[acc]
                            end
                            if acc == CV_COLLISION_ENERGY
                                val = attribute(cv, "value")
                                if val !== nothing
                                    ce = val
                                end
                            end
                        end
                    end
                    if actMethod != ""
                        filter *= "  " * actMethod
                        if ce != ""
                            filter *= "(CE=" * ce * ")"
                        end
                    end
                end
            end
        end

        if !(filter in seen)
            push!(seen, filter)
            push!(info, filter)
        end
    end

    free(xdoc)
    return info
end


"""
    load_mzml_all(filename::String)
Load all spectra from an mzML file. Returns a Vector{MSscan}.
"""
function load_mzml_all(filename::String)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)

    run_elem = find_element(mzml, "run")
    if run_elem === nothing
        free(xdoc)
        error("No <run> element found in mzML file.")
    end

    specList = find_element(run_elem, "spectrumList")
    if specList === nothing
        free(xdoc)
        error("No <spectrumList> element found in mzML file.")
    end

    countStr  = attribute(specList, "count")
    scanCount = countStr !== nothing ? parse(Int, countStr) : 0

    # Buffer with abstract element type — narrowed below. Also remember
    # whether each spectrum carried the MassJ "saved_as_scalar" marker so we
    # can unwrap to a bare value if exactly one was marked.
    buf       = Vector{MScontainer}(undef, scanCount)
    scalar_of = Vector{Bool}(undef, scanCount)
    index = 1
    for spec in child_elements(specList)
        if name(spec) == "spectrum"
            scalar_of[index] = _mzml_is_saved_as_scalar(spec)
            buf[index]       = load_mzml_spectrum(spec, index)
            index += 1
        end
    end

    free(xdoc)
    raw         = buf[1:index-1]
    scalar_vec  = scalar_of[1:index-1]

    # Bit-symmetric scalar round-trip: a file with exactly one spectrum that
    # was saved from a bare MSscan or MSscans returns the bare value.
    if length(raw) == 1 && scalar_vec[1]
        return raw[1]
    end

    # Otherwise narrow Vector{MScontainer} to its concrete element type.
    if isempty(raw)
        return Vector{MSscan}()
    elseif all(x -> x isa MSscans, raw)
        return convert(Vector{MSscans}, raw)
    elseif all(x -> x isa MSscan, raw)
        return convert(Vector{MSscan}, raw)
    else
        return raw  # mixed (rare)
    end
end


"""
    _mzml_is_saved_as_scalar(spec::XMLElement) -> Bool
True when the spectrum carries the MassJ "saved_as_scalar" `userParam`,
written by [`save_mzml`](@ref) for a bare `MSscan` / `MSscans` input.
"""
function _mzml_is_saved_as_scalar(spec::XMLElement)
    for child in child_elements(spec)
        if name(child) == "userParam" &&
           attribute(child, "name") == MASSJ_SCALAR_PARAM &&
           attribute(child, "value") == "true"
            return true
        end
    end
    return false
end


"""
    load_mzml(filename::String, target_num::Int)
Load a single scan from an mzML file by scan number.
"""
function load_mzml(filename::String, target_num::Int)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)

    run_elem = find_element(mzml, "run")
    specList = find_element(run_elem, "spectrumList")

    index = 1
    for spec in child_elements(specList)
        if name(spec) == "spectrum"
            if index == target_num
                scan = load_mzml_spectrum(spec, index)
                free(xdoc)
                return scan
            end
            index += 1
        end
    end

    free(xdoc)
    error("Scan number $target_num not found in mzML file.")
end


"""
    load_mzml_spectrum(spec::XMLElement, scan_index::Int)
Parse a single <spectrum> element and return an MSscan.
"""
function load_mzml_spectrum(spec::XMLElement, scan_index::Int)
    # MS level
    msLevel = parse(Int, get_cv_value(spec, CV_MS_LEVEL, "0"))

    # Spectrum type
    spectrumType = :unknown
    if has_cv_param(spec, CV_CENTROID)
        spectrumType = :centroid
    elseif has_cv_param(spec, CV_PROFILE)
        spectrumType = :profile
    end

    # Polarity
    polarity = ""
    if has_cv_param(spec, CV_POSITIVE_SCAN)
        polarity = "+"
    elseif has_cv_param(spec, CV_NEGATIVE_SCAN)
        polarity = "-"
    end

    # TIC, base peak
    tic = parse(Float64, get_cv_value(spec, CV_TIC, "0.0"))
    basePeakMz = parse(Float64, get_cv_value(spec, CV_BASE_PEAK_MZ, "0.0"))
    basePeakIntensity = parse(Float64, get_cv_value(spec, CV_BASE_PEAK_INT, "0.0"))

    # Retention time from <scanList>/<scan>
    rt = 0.0
    driftTime = -1.0
    compensationVoltage = 0.0
    mobilityType = :none

    scanListElem = find_element(spec, "scanList")
    if scanListElem !== nothing
        for scanElem in child_elements(scanListElem)
            if name(scanElem) == "scan"
                # Retention time
                rtParam = get_cv_param(scanElem, CV_SCAN_START_TIME)
                if rtParam !== nothing
                    rtVal = parse(Float64, attribute(rtParam, "value"))
                    unitAcc = attribute(rtParam, "unitAccession")
                    if unitAcc == CV_UNIT_MINUTE
                        rt = rtVal  # keep in minutes
                    elseif unitAcc == CV_UNIT_SECOND
                        rt = rtVal / 60.0  # convert to minutes
                    else
                        rt = rtVal  # assume minutes
                    end
                end

                # Drift time
                dtParam = get_cv_param(scanElem, CV_DRIFT_TIME)
                if dtParam !== nothing
                    driftTime = parse(Float64, attribute(dtParam, "value"))
                    mobilityType = :DTIMS
                end
                k0Param = get_cv_param(scanElem, CV_INV_K0)
                if k0Param !== nothing
                    driftTime = parse(Float64, attribute(k0Param, "value"))
                    mobilityType = :TIMS
                end

                # Compensation voltage (FAIMS/DMS)
                cvParam = get_cv_param(scanElem, CV_FAIMS_CV)
                if cvParam !== nothing
                    compensationVoltage = parse(Float64, attribute(cvParam, "value"))
                    mobilityType = :FAIMS
                end
                cvParam2 = get_cv_param(scanElem, CV_SELEXION_CV)
                if cvParam2 !== nothing
                    compensationVoltage = parse(Float64, attribute(cvParam2, "value"))
                    mobilityType = :FAIMS
                end

                break  # use first scan element
            end
        end
    end

    # Precursor info
    precursorMz = 0.0
    chargeState = 0
    activationMethod = ""
    collisionEnergy = 0.0

    precList = find_element(spec, "precursorList")
    if precList !== nothing
        for prec in child_elements(precList)
            if name(prec) != "precursor"
                continue
            end

            # Selected ion
            selIonList = find_element(prec, "selectedIonList")
            if selIonList !== nothing
                for selIon in child_elements(selIonList)
                    if name(selIon) == "selectedIon"
                        pmz = get_cv_value(selIon, CV_SELECTED_ION_MZ, "0.0")
                        precursorMz = parse(Float64, pmz)
                        cs = get_cv_value(selIon, CV_CHARGE_STATE, "0")
                        chargeState = parse(Int, cs)
                        break
                    end
                end
            end

            # Activation
            actElem = find_element(prec, "activation")
            if actElem !== nothing
                for cv in child_elements(actElem)
                    if name(cv) == "cvParam"
                        acc = attribute(cv, "accession")
                        if acc !== nothing && haskey(ACTIVATION_METHODS, acc)
                            activationMethod = ACTIVATION_METHODS[acc]
                        end
                        if acc == CV_COLLISION_ENERGY
                            val = attribute(cv, "value")
                            if val !== nothing
                                collisionEnergy = parse(Float64, val)
                            end
                        end
                    end
                end
            end

            break  # use first precursor
        end
    end

    # Binary data arrays
    mz = Float64[]
    int_arr = Float64[]

    bdaList = find_element(spec, "binaryDataArrayList")
    if bdaList !== nothing
        for bda in child_elements(bdaList)
            if name(bda) != "binaryDataArray"
                continue
            end

            is_mz = has_cv_param(bda, CV_MZ_ARRAY)
            is_int = has_cv_param(bda, CV_INT_ARRAY)

            if !is_mz && !is_int
                continue  # skip other arrays (e.g., charge array)
            end

            # Determine precision
            is_64bit = has_cv_param(bda, CV_64BIT)

            # Determine compression
            is_zlib = has_cv_param(bda, CV_ZLIB)

            # Get binary content
            binElem = find_element(bda, "binary")
            if binElem === nothing
                continue
            end
            binContent = content(binElem)
            if isempty(strip(binContent))
                continue
            end

            # Decode: base64 -> optional zlib -> reinterpret
            data = decode(Base64, binContent)
            if is_zlib
                data = Libz.inflate(data)
            end

            if is_64bit
                arr = reinterpret(Float64, data)
            else
                arr = reinterpret(Float32, data)
            end

            # mzML uses little-endian; convert from little-endian to host
            arr = ltoh.(arr)
            arr = convert(Vector{Float64}, arr)

            if is_mz
                mz = arr
            elseif is_int
                int_arr = arr
            end
        end
    end

    scan = MSscan(scan_index, rt, tic, mz, int_arr, msLevel,
                  basePeakMz, basePeakIntensity, precursorMz, polarity,
                  activationMethod, collisionEnergy,
                  chargeState, spectrumType, driftTime, compensationVoltage,
                  mobilityType, Dict{String,Any}())

    # Promote to MSscans if the MassJ export marker is present.
    if _mzml_is_msscans(spec)
        variance = _mzml_extract_variance(spec)
        return _msscans_from_userParams(spec, scan, variance)
    end
    return scan
end


"""
    _msscans_from_userParams(spec::XMLElement, scan::MSscan,
                             variance::Vector{Float64}) -> MSscans
Build an `MSscans` from the vector-valued provenance fields stored as MassJ
`userParam` children of `spec`. Any field missing from the file falls back to
the scalar from `scan` wrapped in a 1-element vector (so partially-formed
files still load cleanly).
"""
function _msscans_from_userParams(spec::XMLElement, scan::MSscan,
                                  variance::Vector{Float64})
    num     = _get_vec_userParam(spec, "MassJ:num",                 Int)
    rtvec   = _get_vec_userParam(spec, "MassJ:rt",                  Float64)
    level   = _get_vec_userParam(spec, "MassJ:level",               Int)
    precvec = _get_vec_userParam(spec, "MassJ:precursor",           Float64)
    pol     = _get_vec_userParam(spec, "MassJ:polarity",            String)
    am      = _get_vec_userParam(spec, "MassJ:activationMethod",    String)
    ce      = _get_vec_userParam(spec, "MassJ:collisionEnergy",     Float64)
    chg     = _get_vec_userParam(spec, "MassJ:chargeState",         Int)
    dt      = _get_vec_userParam(spec, "MassJ:driftTime",           Float64)
    cv      = _get_vec_userParam(spec, "MassJ:compensationVoltage", Float64)

    return MSscans(
        something(num,     [scan.num]),
        something(rtvec,   [scan.rt]),
        scan.tic, scan.mz, scan.int,
        something(level,   [scan.level]),
        scan.basePeakMz, scan.basePeakIntensity,
        something(precvec, [scan.precursor]),
        something(pol,     [scan.polarity]),
        something(am,      [scan.activationMethod]),
        something(ce,      [scan.collisionEnergy]),
        variance,
        something(chg,     [scan.chargeState]),
        scan.spectrumType,
        something(dt,      [scan.driftTime]),
        something(cv,      [scan.compensationVoltage]),
        scan.mobilityType,
        scan.metadata,
    )
end


"""
    _get_vec_userParam(spec, name, T) -> Vector{T} or `nothing`
Read a `userParam` named `name` from `spec` and parse its value as a
pipe-separated vector of `T`. Returns `nothing` if the `userParam` is absent.
"""
function _get_vec_userParam(spec::XMLElement, pname::String, ::Type{T}) where T
    for child in child_elements(spec)
        name(child) == "userParam" || continue
        attribute(child, "name") == pname || continue
        val = attribute(child, "value")
        (val === nothing || isempty(val)) && return T[]
        parts = split(val, '|')
        return T === String ? String.(parts) : parse.(T, parts)
    end
    return nothing
end


"""
    _mzml_is_msscans(spec::XMLElement) -> Bool
True when the `<spectrum>` carries the MassJ "container_type = MSscans"
`userParam` marker (written by [`save_mzml`](@ref) for averaged spectra).
"""
function _mzml_is_msscans(spec::XMLElement)
    for child in child_elements(spec)
        if name(child) == "userParam" &&
           attribute(child, "name") == MASSJ_CONTAINER_PARAM &&
           attribute(child, "value") == "MSscans"
            return true
        end
    end
    return false
end


"""
    _mzml_extract_variance(spec::XMLElement) -> Vector{Float64}
Return the per-m/z variance array stored alongside m/z and intensity by
[`save_mzml`](@ref) for an averaged spectrum. Looks for a `<binaryDataArray>`
carrying the MassJ `userParam name = "MassJ:variance_array"` marker. Returns
an empty vector if the array is absent.
"""
function _mzml_extract_variance(spec::XMLElement)
    bdaList = find_element(spec, "binaryDataArrayList")
    bdaList === nothing && return Float64[]
    for bda in child_elements(bdaList)
        name(bda) == "binaryDataArray" || continue
        is_variance = false
        for child in child_elements(bda)
            if name(child) == "userParam" &&
               attribute(child, "name") == MASSJ_VARIANCE_PARAM
                is_variance = true
                break
            end
        end
        is_variance || continue

        is_64bit = has_cv_param(bda, CV_64BIT)
        is_zlib  = has_cv_param(bda, CV_ZLIB)
        binElem  = find_element(bda, "binary")
        binElem === nothing && return Float64[]
        binContent = content(binElem)
        isempty(strip(binContent)) && return Float64[]

        data = decode(Base64, binContent)
        if is_zlib
            data = Libz.inflate(data)
        end
        arr = is_64bit ? reinterpret(Float64, data) : reinterpret(Float32, data)
        arr = ltoh.(arr)
        return convert(Vector{Float64}, arr)
    end
    return Float64[]
end


"""
    _promote_to_msscans(scan::MSscan, variance::Vector{Float64}) -> MSscans
Wrap a single [`MSscan`](@ref) into an [`MSscans`](@ref) by promoting the
scalar provenance fields to length-1 vectors and attaching `variance` as the
Welford `s` accumulator. Used by [`load`](@ref) when it sees a MassJ MSscans
marker in an mzML / mzXML file.
"""
function _promote_to_msscans(scan::MSscan, variance::Vector{Float64})
    return MSscans(
        [scan.num],
        [scan.rt],
        scan.tic, scan.mz, scan.int,
        [scan.level],
        scan.basePeakMz, scan.basePeakIntensity,
        [scan.precursor],
        [scan.polarity],
        [scan.activationMethod],
        [scan.collisionEnergy],
        variance,
        [scan.chargeState],
        scan.spectrumType,
        [scan.driftTime],
        [scan.compensationVoltage],
        scan.mobilityType,
        scan.metadata,
    )
end


"""
    retention_time_mzml(filename::String)
Extract retention times from an mzML file without loading full spectra.
"""
function retention_time_mzml(filename::String)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)
    run_elem = find_element(mzml, "run")
    specList = find_element(run_elem, "spectrumList")

    rt = Vector{Float64}(undef, 0)

    for spec in child_elements(specList)
        if name(spec) != "spectrum"
            continue
        end
        scanListElem = find_element(spec, "scanList")
        if scanListElem === nothing
            continue
        end
        for scanElem in child_elements(scanListElem)
            if name(scanElem) == "scan"
                rtParam = get_cv_param(scanElem, CV_SCAN_START_TIME)
                if rtParam !== nothing
                    rtVal = parse(Float64, attribute(rtParam, "value"))
                    unitAcc = attribute(rtParam, "unitAccession")
                    if unitAcc == CV_UNIT_SECOND
                        rtVal /= 60.0
                    end
                    push!(rt, rtVal)
                end
                break
            end
        end
    end

    free(xdoc)
    return rt
end
