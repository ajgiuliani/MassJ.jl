"""
Interface to the imzML file format for imaging mass spectrometry.
imzML consists of an XML metadata file (.imzML) and a separate binary data file (.ibd).
The XML follows the mzML schema with additional IMS CV terms for spatial coordinates
and external binary references.
"""


# IMS CV accession constants
const CV_IMS_CONTINUOUS    = "IMS:1000030"
const CV_IMS_PROCESSED     = "IMS:1000031"
const CV_IMS_POSITION_X    = "IMS:1000050"
const CV_IMS_POSITION_Y    = "IMS:1000051"
const CV_IMS_POSITION_Z    = "IMS:1000052"
const CV_IMS_EXTERNAL_DATA = "IMS:1000101"
const CV_IMS_EXT_OFFSET    = "IMS:1000102"
const CV_IMS_EXT_LENGTH    = "IMS:1000103"
const CV_IMS_EXT_ENC_LEN   = "IMS:1000104"
const CV_IMS_UUID          = "IMS:1000080"
const CV_IMS_PIXEL_SIZE_X  = "IMS:1000046"
const CV_IMS_PIXEL_SIZE_Y  = "IMS:1000047"
const CV_IMS_MAX_PIXELS_X  = "IMS:1000042"
const CV_IMS_MAX_PIXELS_Y  = "IMS:1000043"


"""
    parse_referenceable_param_groups(mzml::XMLElement)
Parse the `<referenceableParamGroupList>` element and return a dictionary
mapping group IDs to vectors of `(accession, value)` tuples.
Real-world imzML files (e.g. ProteoWizard output) commonly use these groups
to avoid repeating cvParams on every `<binaryDataArray>`.
"""
function parse_referenceable_param_groups(mzml::XMLElement)
    groups = Dict{String, Vector{Tuple{String,String}}}()
    rpgList = find_element(mzml, "referenceableParamGroupList")
    if rpgList === nothing
        return groups
    end
    for rpg in child_elements(rpgList)
        if name(rpg) != "referenceableParamGroup"
            continue
        end
        gid = attribute(rpg, "id")
        if gid === nothing
            continue
        end
        params = Tuple{String,String}[]
        for cv in child_elements(rpg)
            if name(cv) == "cvParam"
                acc = attribute(cv, "accession")
                val = attribute(cv, "value")
                if acc !== nothing
                    push!(params, (acc, val === nothing ? "" : val))
                end
            end
        end
        groups[gid] = params
    end
    return groups
end


"""
    has_cv_param_resolved(elem::XMLElement, accession::String,
                          ref_groups::Dict{String, Vector{Tuple{String,String}}})
Check whether `elem` has a cvParam with the given accession, either as a direct
child or via a `<referenceableParamGroupRef>`.
"""
function has_cv_param_resolved(elem::XMLElement, accession::String,
                               ref_groups::Dict{String, Vector{Tuple{String,String}}})
    # Check direct cvParam children first
    if has_cv_param(elem, accession)
        return true
    end
    # Check referenced param groups
    for child in child_elements(elem)
        if name(child) == "referenceableParamGroupRef"
            ref = attribute(child, "ref")
            if ref !== nothing && haskey(ref_groups, ref)
                for (acc, _) in ref_groups[ref]
                    if acc == accession
                        return true
                    end
                end
            end
        end
    end
    return false
end


"""
    load_imzml_all(filename::String)
Load all spectra from an imzML file. Returns a `Vector{MSscan}`.
The spatial coordinates (x, y) are stored in the `metadata` dict of each scan.
Requires the companion `.ibd` file in the same directory.
"""
function load_imzml_all(filename::String)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)

    # Parse referenceableParamGroups (used by ProteoWizard-generated imzML)
    ref_groups = parse_referenceable_param_groups(mzml)

    # Determine .ibd file path
    basepath = filename[1:findlast('.', filename)-1]
    ibd_path = basepath * ".ibd"
    if !isfile(ibd_path)
        # Try case variations
        for ext in (".ibd", ".IBD", ".Ibd")
            candidate = basepath * ext
            if isfile(candidate)
                ibd_path = candidate
                break
            end
        end
    end
    if !isfile(ibd_path)
        free(xdoc)
        error("Companion .ibd file not found for $filename")
    end

    # Determine storage mode (continuous vs processed)
    is_continuous = false
    fileDesc = find_element(mzml, "fileDescription")
    if fileDesc !== nothing
        fileContent = find_element(fileDesc, "fileContent")
        if fileContent !== nothing
            if has_cv_param(fileContent, CV_IMS_CONTINUOUS)
                is_continuous = true
            end
        end
    end

    # Get spectrum list
    run_elem = find_element(mzml, "run")
    if run_elem === nothing
        free(xdoc)
        error("No <run> element found in imzML file.")
    end

    specList = find_element(run_elem, "spectrumList")
    if specList === nothing
        free(xdoc)
        error("No <spectrumList> element found in imzML file.")
    end

    countStr = attribute(specList, "count")
    scanCount = countStr !== nothing ? parse(Int, countStr) : 0

    scans = Vector{MSscan}(undef, scanCount)

    # Open .ibd and parse spectra
    open(ibd_path, "r") do ibd_io
        # Skip the 16-byte UUID at the start
        skip(ibd_io, 16)

        # In continuous mode, read the shared m/z array on first encounter
        shared_mz = Float64[]
        shared_mz_read = false

        index = 1
        for spec in child_elements(specList)
            if name(spec) != "spectrum"
                continue
            end

            scan = load_imzml_spectrum(spec, index, ibd_io,
                                       is_continuous, shared_mz, shared_mz_read,
                                       ref_groups)
            scans[index] = scan

            # After first spectrum in continuous mode, cache the shared m/z
            if is_continuous && !shared_mz_read && !isempty(scan.mz)
                append!(shared_mz, scan.mz)
                shared_mz_read = true
            end

            index += 1
        end
    end

    free(xdoc)
    return scans[1:min(scanCount, length(scans))]
end


"""
    load_imzml_spectrum(spec::XMLElement, scan_index::Int, ibd_io::IO,
                        is_continuous::Bool, shared_mz::Vector{Float64},
                        shared_mz_read::Bool,
                        ref_groups::Dict{String, Vector{Tuple{String,String}}})
Parse a single <spectrum> element from imzML and read its binary data from the .ibd file.
Resolves `referenceableParamGroupRef` elements via `ref_groups`.
"""
function load_imzml_spectrum(spec::XMLElement, scan_index::Int, ibd_io::IO,
                              is_continuous::Bool, shared_mz::Vector{Float64},
                              shared_mz_read::Bool,
                              ref_groups::Dict{String, Vector{Tuple{String,String}}}=Dict{String, Vector{Tuple{String,String}}}())
    # MS level
    msLevel = parse(Int, get_cv_value(spec, CV_MS_LEVEL, "1"))

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

    # Spatial coordinates and retention time
    rt = 0.0
    pos_x = 0
    pos_y = 0
    pos_z = 0

    scanListElem = find_element(spec, "scanList")
    if scanListElem !== nothing
        for scanElem in child_elements(scanListElem)
            if name(scanElem) == "scan"
                # Retention time
                rtParam = get_cv_param(scanElem, CV_SCAN_START_TIME)
                if rtParam !== nothing
                    rtVal = parse(Float64, attribute(rtParam, "value"))
                    unitAcc = attribute(rtParam, "unitAccession")
                    if unitAcc == CV_UNIT_SECOND
                        rtVal /= 60.0
                    end
                    rt = rtVal
                end

                # Spatial coordinates
                xParam = get_cv_param(scanElem, CV_IMS_POSITION_X)
                if xParam !== nothing
                    pos_x = parse(Int, attribute(xParam, "value"))
                end
                yParam = get_cv_param(scanElem, CV_IMS_POSITION_Y)
                if yParam !== nothing
                    pos_y = parse(Int, attribute(yParam, "value"))
                end
                zParam = get_cv_param(scanElem, CV_IMS_POSITION_Z)
                if zParam !== nothing
                    pos_z = parse(Int, attribute(zParam, "value"))
                end

                break
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
            break
        end
    end

    # Read binary data arrays from .ibd file
    mz = Float64[]
    int_arr = Float64[]

    bdaList = find_element(spec, "binaryDataArrayList")
    if bdaList !== nothing
        for bda in child_elements(bdaList)
            if name(bda) != "binaryDataArray"
                continue
            end

            # Check array type using both direct cvParams and referenced groups
            is_mz = has_cv_param_resolved(bda, CV_MZ_ARRAY, ref_groups)
            is_int = has_cv_param_resolved(bda, CV_INT_ARRAY, ref_groups)
            if !is_mz && !is_int
                continue
            end

            # Check if external data
            if !has_cv_param_resolved(bda, CV_IMS_EXTERNAL_DATA, ref_groups)
                continue
            end

            # Get offset and length
            offset_str = get_cv_value(bda, CV_IMS_EXT_OFFSET, "0")
            arr_len_str = get_cv_value(bda, CV_IMS_EXT_LENGTH, "0")
            enc_len_str = get_cv_value(bda, CV_IMS_EXT_ENC_LEN, "0")

            offset = parse(Int, offset_str)
            arr_len = parse(Int, arr_len_str)
            enc_len = parse(Int, enc_len_str)

            if arr_len == 0
                continue
            end

            # In continuous mode, reuse shared m/z array
            if is_continuous && is_mz && shared_mz_read
                mz = copy(shared_mz)
                continue
            end

            # Determine precision (check both direct and referenced)
            is_64bit = has_cv_param_resolved(bda, CV_64BIT, ref_groups)

            # Read from .ibd file
            seek(ibd_io, offset)
            raw_data = read(ibd_io, enc_len)

            # Check for zlib compression (check both direct and referenced)
            if has_cv_param_resolved(bda, CV_ZLIB, ref_groups)
                raw_data = Libz.inflate(raw_data)
            end

            if is_64bit
                arr = reinterpret(Float64, raw_data)
            else
                arr = reinterpret(Float32, raw_data)
            end

            # imzML uses little-endian (same as host on x86)
            arr = ltoh.(arr)
            arr = convert(Vector{Float64}, copy(arr))

            if is_mz
                mz = arr
            elseif is_int
                int_arr = arr
            end
        end
    end

    # Compute TIC and base peak if not in XML
    if tic == 0.0 && !isempty(int_arr)
        tic = sum(int_arr)
    end
    if basePeakIntensity == 0.0 && !isempty(int_arr)
        maxidx = argmax(int_arr)
        basePeakMz = mz[maxidx]
        basePeakIntensity = int_arr[maxidx]
    end

    # Store spatial coordinates in metadata
    metadata = Dict{String,Any}(
        "position_x" => pos_x,
        "position_y" => pos_y,
    )
    if pos_z != 0
        metadata["position_z"] = pos_z
    end

    return MSscan(scan_index, rt, tic, mz, int_arr, msLevel,
                  basePeakMz, basePeakIntensity, precursorMz, polarity,
                  activationMethod, collisionEnergy,
                  chargeState, spectrumType, -1.0, 0.0, :none, metadata)
end


"""
    info_imzml(filename::String, info::Vector{String}, verbose::Bool=false)
Returns summary information about an imzML file.
"""
function info_imzml(filename::String, info::Vector{String}, verbose::Bool=false)
    xdoc = parse_file(filename)
    mzml = find_mzml_root(xdoc)

    # Storage mode
    is_continuous = false
    fileDesc = find_element(mzml, "fileDescription")
    if fileDesc !== nothing
        fileContent = find_element(fileDesc, "fileContent")
        if fileContent !== nothing
            if has_cv_param(fileContent, CV_IMS_CONTINUOUS)
                is_continuous = true
            end
        end
    end

    if verbose
        push!(info, "Storage mode: " * (is_continuous ? "continuous" : "processed"))

        # Instrument info (same as mzML)
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
                end
            end
        end
    end

    # Count spectra and collect scan info
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
        push!(info, scanCount * " spectra")
    end

    # Determine image dimensions from scan settings or from scanning all spectra
    max_x = 0
    max_y = 0
    seen = Set{String}()

    for spec in child_elements(specList)
        if name(spec) != "spectrum"
            continue
        end

        # Track spatial extent
        scanListElem = find_element(spec, "scanList")
        if scanListElem !== nothing
            for scanElem in child_elements(scanListElem)
                if name(scanElem) == "scan"
                    xParam = get_cv_param(scanElem, CV_IMS_POSITION_X)
                    if xParam !== nothing
                        x = parse(Int, attribute(xParam, "value"))
                        max_x = max(max_x, x)
                    end
                    yParam = get_cv_param(scanElem, CV_IMS_POSITION_Y)
                    if yParam !== nothing
                        y = parse(Int, attribute(yParam, "value"))
                        max_y = max(max_y, y)
                    end
                    break
                end
            end
        end

        # Collect unique scan descriptions
        msLevel = get_cv_value(spec, CV_MS_LEVEL, "1")
        desc = "MS" * msLevel
        if has_cv_param(spec, CV_POSITIVE_SCAN)
            desc *= "+"
        elseif has_cv_param(spec, CV_NEGATIVE_SCAN)
            desc *= "-"
        end
        if !(desc in seen)
            push!(seen, desc)
            push!(info, desc)
        end
    end

    if max_x > 0 || max_y > 0
        push!(info, "Image dimensions: " * string(max_x) * " x " * string(max_y) * " pixels")
    end

    free(xdoc)
    return info
end
