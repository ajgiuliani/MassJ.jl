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
const CV_TIME_ARRAY        = "MS:1000595"
const CV_TIC_CHROMATOGRAM  = "MS:1000235"
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

# Optional cvParams populated into MSscan.metadata when present in the source
# mzML. Locations: see _read_mzml_extra_metadata.
const CV_SPECTRUM_TITLE       = "MS:1000796"   # in <spectrum>
const CV_LOWEST_OBSERVED_MZ   = "MS:1000528"   # in <spectrum>
const CV_HIGHEST_OBSERVED_MZ  = "MS:1000527"   # in <spectrum>
const CV_MASS_RESOLVING_POWER = "MS:1000800"   # in <scan>
const CV_ION_INJECTION_TIME   = "MS:1000927"   # in <scan>
const CV_SCAN_WINDOW_LOWER    = "MS:1000501"   # in <scanWindow>
const CV_SCAN_WINDOW_UPPER    = "MS:1000500"   # in <scanWindow>
const CV_UNIT_MILLISECOND     = "UO:0000028"
const CV_UNIT_MZ              = "MS:1000040"
const CV_UNIT_DA              = "UO:0000221"

# Per-spectrum cvParams (Phase 1 — feature/mzml-full-roundtrip)
const CV_MSN_SPECTRUM         = "MS:1000580"   # <spectrum> marker (level>=2)
const CV_FILTER_STRING        = "MS:1000512"   # in <scan>, Thermo-specific
const CV_SELECTED_PEAK_INT    = "MS:1000042"   # in <selectedIon>
const CV_ISO_WINDOW_TARGET    = "MS:1000827"   # in <isolationWindow>
const CV_ISO_WINDOW_LOWER     = "MS:1000828"   # in <isolationWindow>
const CV_ISO_WINDOW_UPPER     = "MS:1000829"   # in <isolationWindow>

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

    # File-level metadata (instrument, software, source file, processing).
    file_md = _read_mzml_file_metadata(mzml)

    run_elem = find_element(mzml, "run")
    if run_elem === nothing
        free(xdoc)
        error("No <run> element found in mzML file.")
    end

    # Pre-computed chromatograms (TIC etc.) — sibling of <spectrumList>.
    chromatograms = _read_mzml_chromatograms(run_elem)

    specList = find_element(run_elem, "spectrumList")
    if specList === nothing
        free(xdoc)
        error("No <spectrumList> element found in mzML file.")
    end

    countStr  = attribute(specList, "count")
    scanCount = countStr !== nothing ? parse(Int, countStr) : 0

    # Run-level default instrument configuration, used when a scan omits its own ref.
    default_ic = attribute(run_elem, "defaultInstrumentConfigurationRef")

    # Remember whether each spectrum carried the MassJ "saved_as_scalar" marker
    # so we can unwrap to a bare value if exactly one was marked.
    buf       = Vector{MSscans}(undef, scanCount)
    scalar_of = Vector{Bool}(undef, scanCount)
    index = 1
    for spec in child_elements(specList)
        if name(spec) == "spectrum"
            scalar_of[index] = _mzml_is_saved_as_scalar(spec)
            buf[index]       = load_mzml_spectrum(spec, index; default_ic_ref = default_ic)
            index += 1
        end
    end

    free(xdoc)
    raw         = buf[1:index-1]
    scalar_vec  = scalar_of[1:index-1]

    # Bit-symmetric scalar round-trip: a file with exactly one spectrum that was
    # saved from a bare spectrum returns the bare value (no MSrun wrapper, no
    # file-level metadata).
    if length(raw) == 1 && scalar_vec[1]
        return raw[1]
    end

    # Multi-spectrum mzML files — wrap in MSrun so the file-level metadata and
    # any pre-computed chromatograms ride along. MSrun is an
    # AbstractVector{MSscans}, so it behaves like the spectrum vector itself.
    return MSrun(raw, file_md, chromatograms)
end


"""
    _mzml_is_saved_as_scalar(spec::XMLElement) -> Bool
True when the spectrum carries the MassJ "saved_as_scalar" `userParam`,
written by [`save_mzml`](@ref) for a bare spectrum input.
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
    default_ic = attribute(run_elem, "defaultInstrumentConfigurationRef")

    index = 1
    for spec in child_elements(specList)
        if name(spec) == "spectrum"
            if index == target_num
                scan = load_mzml_spectrum(spec, index; default_ic_ref = default_ic)
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
    load_mzml_spectrum(spec::XMLElement, scan_index::Int; default_ic_ref = nothing)
Parse a single <spectrum> element and return an MSscans. `default_ic_ref` supplies the run's
`defaultInstrumentConfigurationRef` used when a scan omits its own `instrumentConfigurationRef`.
"""
function load_mzml_spectrum(spec::XMLElement, scan_index::Int;
                            default_ic_ref::Union{AbstractString,Nothing} = nothing)
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
    ic_ref = nothing   # per-scan instrumentConfigurationRef

    scanListElem = find_element(spec, "scanList")
    if scanListElem !== nothing
        for scanElem in child_elements(scanListElem)
            if name(scanElem) == "scan"
                # Instrument configuration reference (falls back to run default below)
                ic_ref = attribute(scanElem, "instrumentConfigurationRef")

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

    extra_md = _read_mzml_extra_metadata(spec)
    ref = ic_ref !== nothing ? ic_ref : default_ic_ref
    ref !== nothing && (extra_md["instrument_configuration_ref"] = String(ref))

    scan = MSscans(scan_index, rt, tic, mz, int_arr, msLevel,
                   basePeakMz, basePeakIntensity, precursorMz, polarity,
                   activationMethod, collisionEnergy,
                   chargeState, spectrumType, driftTime, compensationVoltage,
                   mobilityType, extra_md)

    # Reconstruct a composite spectrum if the MassJ export marker is present.
    if _mzml_is_msscans(spec)
        variance = _mzml_extract_variance(spec)
        return _msscans_from_userParams(spec, scan, variance)
    end
    return scan
end


"""
    _msscans_from_userParams(spec::XMLElement, scan::MSscans,
                             variance::Vector{Float64}) -> MSscans
Rebuild a composite `MSscans` from the vector-valued provenance fields stored as
MassJ `userParam` children of `spec`. Any field missing from the file falls back
to the corresponding provenance vector of the single-scan base `scan` (so
partially-formed files still load cleanly).
"""
function _msscans_from_userParams(spec::XMLElement, scan::MSscans,
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
        something(num,     scan.num),
        something(rtvec,   scan.rt),
        scan.tic, scan.mz, scan.int,
        something(level,   scan.level),
        scan.basePeakMz, scan.basePeakIntensity,
        something(precvec, scan.precursor),
        something(pol,     scan.polarity),
        something(am,      scan.activationMethod),
        something(ce,      scan.collisionEnergy),
        variance,
        something(chg,     scan.chargeState),
        scan.spectrumType,
        something(dt,      scan.driftTime),
        something(cv,      scan.compensationVoltage),
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



"""
    _read_mzml_extra_metadata(spec::XMLElement) -> Dict{String,Any}
Pull optional cvParams off an mzML `<spectrum>` element and return them as a
`Dict`. The dict is populated only for cvParams that are actually present in
the source file; keys for absent fields are not added.

Keys populated (when present in the source mzML):

| Key                    | mzML cvParam                       | XML location               |
|------------------------|------------------------------------|----------------------------|
| `"spectrum_title"`     | MS:1000796 spectrum title          | direct child of `<spectrum>` |
| `"lowest_observed_mz"` | MS:1000528 lowest observed m/z     | direct child of `<spectrum>` |
| `"highest_observed_mz"`| MS:1000527 highest observed m/z    | direct child of `<spectrum>` |
| `"mass_resolving_power"`| MS:1000800 mass resolving power   | `<scan>` inside `<scanList>` |
| `"ion_injection_time"` | MS:1000927 ion injection time (ms) | `<scan>` inside `<scanList>` |
| `"scan_window_lower"`  | MS:1000501 scan window lower limit | `<scanWindow>` inside `<scan>` |
| `"scan_window_upper"`  | MS:1000500 scan window upper limit | `<scanWindow>` inside `<scan>` |

Used by [`load_mzml_spectrum`](@ref) to populate `MSscan.metadata`.
"""
function _read_mzml_extra_metadata(spec::XMLElement)
    md = Dict{String,Any}()

    # -- Direct cvParam children of <spectrum> ------------------------------
    title = get_cv_param(spec, CV_SPECTRUM_TITLE)
    if title !== nothing
        val = attribute(title, "value")
        if val !== nothing && !isempty(val)
            md["spectrum_title"] = val
        end
    end
    lo = get_cv_param(spec, CV_LOWEST_OBSERVED_MZ)
    if lo !== nothing
        val = attribute(lo, "value")
        if val !== nothing && !isempty(val)
            md["lowest_observed_mz"] = parse(Float64, val)
        end
    end
    hi = get_cv_param(spec, CV_HIGHEST_OBSERVED_MZ)
    if hi !== nothing
        val = attribute(hi, "value")
        if val !== nothing && !isempty(val)
            md["highest_observed_mz"] = parse(Float64, val)
        end
    end

    # -- Inside <scanList>/<scan> -------------------------------------------
    scanListElem = find_element(spec, "scanList")
    if scanListElem !== nothing
        for scanElem in child_elements(scanListElem)
            name(scanElem) == "scan" || continue

            rp = get_cv_param(scanElem, CV_MASS_RESOLVING_POWER)
            if rp !== nothing
                val = attribute(rp, "value")
                if val !== nothing && !isempty(val)
                    md["mass_resolving_power"] = parse(Float64, val)
                end
            end
            it = get_cv_param(scanElem, CV_ION_INJECTION_TIME)
            if it !== nothing
                val = attribute(it, "value")
                if val !== nothing && !isempty(val)
                    md["ion_injection_time"] = parse(Float64, val)
                end
            end
            fs = get_cv_param(scanElem, CV_FILTER_STRING)
            if fs !== nothing
                val = attribute(fs, "value")
                if val !== nothing && !isempty(val)
                    md["filter_string"] = val
                end
            end

            # Scan windows (typically one per scan; we take the first)
            swList = find_element(scanElem, "scanWindowList")
            if swList !== nothing
                for sw in child_elements(swList)
                    name(sw) == "scanWindow" || continue
                    lower = get_cv_param(sw, CV_SCAN_WINDOW_LOWER)
                    if lower !== nothing
                        val = attribute(lower, "value")
                        if val !== nothing && !isempty(val)
                            md["scan_window_lower"] = parse(Float64, val)
                        end
                    end
                    upper = get_cv_param(sw, CV_SCAN_WINDOW_UPPER)
                    if upper !== nothing
                        val = attribute(upper, "value")
                        if val !== nothing && !isempty(val)
                            md["scan_window_upper"] = parse(Float64, val)
                        end
                    end
                    break  # first window only
                end
            end

            break  # first scan only
        end
    end

    # -- Inside <precursorList>/<precursor> ---------------------------------
    # Isolation window (target + offsets) and selected ion peak intensity.
    precList = find_element(spec, "precursorList")
    if precList !== nothing
        for prec in child_elements(precList)
            name(prec) == "precursor" || continue

            isoWin = find_element(prec, "isolationWindow")
            if isoWin !== nothing
                for (acc, key) in (
                        (CV_ISO_WINDOW_TARGET, "isolation_window_target_mz"),
                        (CV_ISO_WINDOW_LOWER,  "isolation_window_lower_offset"),
                        (CV_ISO_WINDOW_UPPER,  "isolation_window_upper_offset"),
                    )
                    cv = get_cv_param(isoWin, acc)
                    if cv !== nothing
                        val = attribute(cv, "value")
                        if val !== nothing && !isempty(val)
                            md[key] = parse(Float64, val)
                        end
                    end
                end
            end

            siList = find_element(prec, "selectedIonList")
            if siList !== nothing
                for si in child_elements(siList)
                    name(si) == "selectedIon" || continue
                    pkInt = get_cv_param(si, CV_SELECTED_PEAK_INT)
                    if pkInt !== nothing
                        val = attribute(pkInt, "value")
                        if val !== nothing && !isempty(val)
                            md["selected_ion_peak_intensity"] = parse(Float64, val)
                        end
                    end
                    break  # first selected ion only
                end
            end

            break  # first precursor only
        end
    end

    return md
end


# ============================================================================
# File-level metadata parsing (for MSrun round-trip)
# ============================================================================

"""
    _cvparams_to_dicts(elem::XMLElement) -> Vector{Dict{String,String}}
Collect all `cvParam` children of `elem` as a vector of `Dict` snapshots.
Each Dict has at least `"accession"` and `"name"`; `value`, `unitAccession`,
`unitName`, `unitCvRef` are included when present in the source XML.
"""
function _cvparams_to_dicts(elem::XMLElement)
    out = Dict{String,String}[]
    for child in child_elements(elem)
        name(child) == "cvParam" || continue
        d = Dict{String,String}()
        for attr in ("accession", "name", "value",
                     "unitAccession", "unitName", "unitCvRef")
            v = attribute(child, attr)
            v !== nothing && (d[attr] = v)
        end
        push!(out, d)
    end
    return out
end


# Collect the `ref` attributes of every <referenceableParamGroupRef> child.
function _param_group_refs(elem::XMLElement)
    refs = String[]
    for child in child_elements(elem)
        if name(child) == "referenceableParamGroupRef"
            r = attribute(child, "ref")
            r !== nothing && push!(refs, r)
        end
    end
    return refs
end


"""
    _read_mzml_file_metadata(mzml::XMLElement) -> Dict{String,Any}
Parse the top-level mzML metadata sections (`fileDescription`,
`referenceableParamGroupList`, `softwareList`, `instrumentConfigurationList`,
`dataProcessingList`) into a nested Dict suitable for storage in
`MSrun.metadata`. The structure is preserved well enough for [`save`](@ref)
to re-emit a semantically equivalent mzML.
"""
function _read_mzml_file_metadata(mzml::XMLElement)
    md = Dict{String,Any}()

    # -- <fileDescription>/<fileContent> + <sourceFileList> -----------------
    fd = find_element(mzml, "fileDescription")
    if fd !== nothing
        fc = find_element(fd, "fileContent")
        if fc !== nothing
            md["file_content"] = _cvparams_to_dicts(fc)
        end
        sfList = find_element(fd, "sourceFileList")
        if sfList !== nothing
            sfs = Dict{String,Any}[]
            for sf in child_elements(sfList)
                name(sf) == "sourceFile" || continue
                d = Dict{String,Any}()
                for attr in ("id", "name", "location")
                    v = attribute(sf, attr)
                    v !== nothing && (d[attr] = v)
                end
                d["cv_params"] = _cvparams_to_dicts(sf)
                push!(sfs, d)
            end
            md["source_files"] = sfs
        end
    end

    # -- <referenceableParamGroupList> --------------------------------------
    pgList = find_element(mzml, "referenceableParamGroupList")
    if pgList !== nothing
        groups = Dict{String,Any}[]
        for pg in child_elements(pgList)
            name(pg) == "referenceableParamGroup" || continue
            d = Dict{String,Any}()
            id = attribute(pg, "id"); id !== nothing && (d["id"] = id)
            d["cv_params"] = _cvparams_to_dicts(pg)
            push!(groups, d)
        end
        md["param_groups"] = groups
    end

    # -- <softwareList> -----------------------------------------------------
    swList = find_element(mzml, "softwareList")
    if swList !== nothing
        sws = Dict{String,Any}[]
        for sw in child_elements(swList)
            name(sw) == "software" || continue
            d = Dict{String,Any}()
            for attr in ("id", "version")
                v = attribute(sw, attr)
                v !== nothing && (d[attr] = v)
            end
            d["cv_params"] = _cvparams_to_dicts(sw)
            push!(sws, d)
        end
        md["software"] = sws
    end

    # -- <instrumentConfigurationList> --------------------------------------
    icList = find_element(mzml, "instrumentConfigurationList")
    if icList !== nothing
        instrs = Dict{String,Any}[]
        for ic in child_elements(icList)
            name(ic) == "instrumentConfiguration" || continue
            d = Dict{String,Any}()
            id = attribute(ic, "id"); id !== nothing && (d["id"] = id)
            d["param_group_refs"] = _param_group_refs(ic)
            d["cv_params"]        = _cvparams_to_dicts(ic)

            comps = Dict{String,Any}[]
            compList = find_element(ic, "componentList")
            if compList !== nothing
                for c in child_elements(compList)
                    cname = name(c)
                    cname ∈ ("source", "analyzer", "detector") || continue
                    cd = Dict{String,Any}("type" => cname)
                    ord = attribute(c, "order")
                    ord !== nothing && (cd["order"] = ord)
                    cd["cv_params"] = _cvparams_to_dicts(c)
                    push!(comps, cd)
                end
            end
            d["components"] = comps

            swRef = find_element(ic, "softwareRef")
            if swRef !== nothing
                r = attribute(swRef, "ref")
                r !== nothing && (d["software_ref"] = r)
            end
            push!(instrs, d)
        end
        md["instruments"] = instrs
    end

    # -- <dataProcessingList> -----------------------------------------------
    dpList = find_element(mzml, "dataProcessingList")
    if dpList !== nothing
        dps = Dict{String,Any}[]
        for dp in child_elements(dpList)
            name(dp) == "dataProcessing" || continue
            d = Dict{String,Any}()
            id = attribute(dp, "id"); id !== nothing && (d["id"] = id)
            methods = Dict{String,Any}[]
            for pm in child_elements(dp)
                name(pm) == "processingMethod" || continue
                m = Dict{String,Any}()
                for attr in ("order", "softwareRef")
                    v = attribute(pm, attr)
                    v !== nothing && (m[attr] = v)
                end
                m["cv_params"] = _cvparams_to_dicts(pm)
                push!(methods, m)
            end
            d["methods"] = methods
            push!(dps, d)
        end
        md["data_processing"] = dps
    end

    return md
end


"""
    _read_mzml_chromatograms(run_elem::XMLElement) -> Vector{IonCurrent}
Parse the `<chromatogramList>` (sibling of `<spectrumList>`) into a vector of
[`IonCurrent`](@ref) traces (axis `:rt`). Only chromatograms carrying the
`MS:1000235` (total ion current chromatogram) marker are loaded; other
chromatogram types (extracted ion, base peak, ...) are skipped with a warning.
Returns an empty vector when the element is absent.

The binary-array decoding pipeline is the same as for spectrum arrays:
base64 → optional zlib → reinterpret as Float64/Float32 → host order.
"""
function _read_mzml_chromatograms(run_elem::XMLElement)
    out = IonCurrent[]
    chList = find_element(run_elem, "chromatogramList")
    chList === nothing && return out

    for chrom in child_elements(chList)
        name(chrom) == "chromatogram" || continue

        if !has_cv_param(chrom, CV_TIC_CHROMATOGRAM)
            @debug "skipping non-TIC chromatogram (MassJ only round-trips TIC chromatograms)" id=attribute(chrom, "id")
            continue
        end

        time_arr  = Float64[]
        int_arr   = Float64[]
        bdaList = find_element(chrom, "binaryDataArrayList")
        bdaList === nothing && continue

        for bda in child_elements(bdaList)
            name(bda) == "binaryDataArray" || continue

            is_time = has_cv_param(bda, CV_TIME_ARRAY)
            is_int  = has_cv_param(bda, CV_INT_ARRAY)
            (is_time || is_int) || continue

            is_64bit = has_cv_param(bda, CV_64BIT)
            is_zlib  = has_cv_param(bda, CV_ZLIB)

            binElem = find_element(bda, "binary")
            binElem === nothing && continue
            binContent = content(binElem)
            isempty(strip(binContent)) && continue

            data = decode(Base64, binContent)
            if is_zlib
                data = Libz.inflate(data)
            end
            arr = is_64bit ? reinterpret(Float64, data) : reinterpret(Float32, data)
            arr = ltoh.(arr)
            arr = convert(Vector{Float64}, arr)

            if is_time
                time_arr = arr
            else
                int_arr = arr
            end
        end

        if !isempty(time_arr) && !isempty(int_arr)
            push!(out, IonCurrent(time_arr, int_arr; axis = :rt))
        end
    end
    return out
end
