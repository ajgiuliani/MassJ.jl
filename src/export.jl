"""
Export MassJ data to standard mass-spectrometry file formats.

The top-level entry point [`save`](@ref) dispatches on the file extension:

* `.mzML`  → [`save_mzml`](@ref)
* `.mzXML` → [`save_mzxml`](@ref)

Both writers round-trip through MassJ's own readers — i.e. `load(save_path)`
recovers the same `MSscan` data. Format-level metadata (`fileDescription`,
`instrumentConfiguration`, `dataProcessing`, the `cvList`) is emitted in a
minimal-but-valid form; downstream tools that need rich provenance fields
should expect them blank.
"""

# User Interface.
# ---------------

export save


"""
    save(data, filename::AbstractString; kwargs...)
Export MassJ data to a file. The format is selected from the extension:

| Extension | Writer            | Notes                                 |
|-----------|-------------------|---------------------------------------|
| `.mzML`   | [`save_mzml`](@ref)  | PSI standard; little-endian arrays    |
| `.mzXML`  | [`save_mzxml`](@ref) | Legacy format; big-endian arrays      |

Accepts a single [`MSscan`](@ref) / [`MSscans`](@ref) or a `Vector{MSscan}`.
Common keyword arguments:

* `precision::Int = 64`   — `64` for `Float64` arrays, `32` for `Float32`
* `compress::Bool = true` — zlib-compress the binary arrays

# Examples
```julia
scans = load("input.mzML")
save(scans, "output.mzML")                     # round-trip
save(scans, "output.mzXML"; precision = 32)    # smaller file, less precision
save(scans[1], "single_scan.mzML")             # one scan
```
"""
function save(data, filename::AbstractString; kwargs...)
    ext = lowercase(splitext(filename)[2])
    ext = startswith(ext, ".") ? ext[2:end] : ext
    if ext == "mzml"
        return save_mzml(filename, data; kwargs...)
    elseif ext == "mzxml"
        return save_mzxml(filename, data; kwargs...)
    else
        error("save: unsupported file format '.$ext' (supported: .mzML, .mzXML)")
    end
end


# Binary encoding shared by mzML / mzXML writers.
# ---------------------------------------------------------------------------

"""
    _encode_binary(data; precision = 64, compress = true, endian = :little)
        -> (base64_string, byte_length)
Encode a real-valued vector to the base64-of-zlib-of-raw-bytes payload expected
by mzML/mzXML `<binary>` / `<peaks>` elements.

`endian` is `:little` for mzML or `:big` for mzXML (network byte order).
"""
function _encode_binary(data::AbstractVector{<:Real};
                        precision::Int  = 64,
                        compress::Bool  = true,
                        endian::Symbol  = :little)
    precision ∈ (32, 64) ||
        error("_encode_binary: precision must be 32 or 64 (got $precision)")
    endian ∈ (:little, :big) ||
        error("_encode_binary: endian must be :little or :big (got :$endian)")

    arr = precision == 64 ? Vector{Float64}(data) : Vector{Float32}(data)
    arr = endian === :little ? htol.(arr) : hton.(arr)
    raw = collect(reinterpret(UInt8, arr))
    bytes = compress ? read(Libz.ZlibDeflateInputStream(raw)) : raw
    b64 = String(Codecs.encode(Codecs.Base64, bytes))
    return b64, length(bytes)
end


# Helper to add a cvParam child element with the usual attributes.
function _cvParam(parent::XMLElement, accession::String, name::String;
                  value::AbstractString = "",
                  unit_cv::AbstractString = "",
                  unit_acc::AbstractString = "",
                  unit_name::AbstractString = "")
    cv = new_child(parent, "cvParam")
    set_attribute(cv, "cvRef", "MS")
    set_attribute(cv, "accession", accession)
    set_attribute(cv, "name", name)
    if !isempty(value)
        set_attribute(cv, "value", value)
    end
    if !isempty(unit_acc)
        set_attribute(cv, "unitCvRef",     unit_cv)
        set_attribute(cv, "unitAccession", unit_acc)
        set_attribute(cv, "unitName",      unit_name)
    end
    return cv
end


# ============================================================================
# mzML writer
# ============================================================================

"""
    save_mzml(filename::AbstractString, data;
              precision::Int = 64, compress::Bool = true) -> filename
Write a [`MSscan`](@ref), [`MSscans`](@ref), or `Vector{MSscan}` to an mzML
file. The emitted file is minimal-but-valid and round-trips through
[`load`](@ref) — m/z, intensity, MS level, polarity, retention time, precursor
m/z, charge state, activation method, and collision energy are preserved.

Optional keywords:
* `precision = 64` — `64` for `Float64` arrays, `32` for `Float32` (smaller, lossy)
* `compress = true` — zlib-compress the binary arrays
"""
function save_mzml(filename::AbstractString, scans::Vector{MSscan};
                   precision::Int = 64, compress::Bool = true)
    return _save_mzml_vector(filename, scans;
                             precision = precision, compress = compress,
                             scalar = false)
end

function save_mzml(filename::AbstractString, scan::MSscan;
                   precision::Int = 64, compress::Bool = true)
    # `scalar = true` records that the input was a bare MSscan, so `load`
    # returns a bare MSscan (not a 1-element Vector) on round-trip.
    return _save_mzml_vector(filename, [scan];
                             precision = precision, compress = compress,
                             scalar = true)
end

function _save_mzml_vector(filename::AbstractString, scans::Vector{MSscan};
                           precision::Int, compress::Bool, scalar::Bool)
    xdoc  = XMLDocument()
    xroot = create_root(xdoc, "mzML")
    set_attribute(xroot, "xmlns",     "http://psi.hupo.org/ms/mzml")
    set_attribute(xroot, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    set_attribute(xroot, "version",   "1.1.0")

    _mzml_cvList(xroot)
    _mzml_fileDescription(xroot)
    _mzml_softwareList(xroot)
    _mzml_instrumentConfigurationList(xroot)
    _mzml_dataProcessingList(xroot)

    run = new_child(xroot, "run")
    set_attribute(run, "id", "run1")
    set_attribute(run, "defaultInstrumentConfigurationRef", "IC1")

    specList = new_child(run, "spectrumList")
    set_attribute(specList, "count", string(length(scans)))
    set_attribute(specList, "defaultDataProcessingRef", "MassJExport")

    for (i, scan) in enumerate(scans)
        _mzml_spectrum(specList, scan, i - 1;
                       precision = precision, compress = compress,
                       scalar = scalar)
    end

    save_file(xdoc, filename)
    free(xdoc)
    return filename
end

function save_mzml(filename::AbstractString, scan::MSscans;
                   precision::Int = 64, compress::Bool = true)
    return _save_mzml_msscans_vector(filename, [scan];
                                     precision = precision, compress = compress,
                                     scalar = true)
end

function save_mzml(filename::AbstractString, scans::Vector{MSscans};
                   precision::Int = 64, compress::Bool = true)
    return _save_mzml_msscans_vector(filename, scans;
                                     precision = precision, compress = compress,
                                     scalar = false)
end

function _save_mzml_msscans_vector(filename::AbstractString, scans::Vector{MSscans};
                                   precision::Int, compress::Bool, scalar::Bool)
    xdoc  = XMLDocument()
    xroot = create_root(xdoc, "mzML")
    set_attribute(xroot, "xmlns",     "http://psi.hupo.org/ms/mzml")
    set_attribute(xroot, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    set_attribute(xroot, "version",   "1.1.0")

    _mzml_cvList(xroot)
    _mzml_fileDescription(xroot)
    _mzml_softwareList(xroot)
    _mzml_instrumentConfigurationList(xroot)
    _mzml_dataProcessingList(xroot)

    run = new_child(xroot, "run")
    set_attribute(run, "id", "run1")
    set_attribute(run, "defaultInstrumentConfigurationRef", "IC1")

    specList = new_child(run, "spectrumList")
    set_attribute(specList, "count", string(length(scans)))
    set_attribute(specList, "defaultDataProcessingRef", "MassJExport")

    for (i, sc) in enumerate(scans)
        _mzml_msscans_spectrum(specList, sc, i - 1;
                               precision = precision, compress = compress,
                               scalar = scalar)
    end

    save_file(xdoc, filename)
    free(xdoc)
    return filename
end


# -- mzML header helpers ------------------------------------------------------

function _mzml_cvList(root::XMLElement)
    cvList = new_child(root, "cvList")
    set_attribute(cvList, "count", "2")
    cv1 = new_child(cvList, "cv")
    set_attribute(cv1, "id",       "MS")
    set_attribute(cv1, "fullName", "Proteomics Standards Initiative Mass Spectrometry Ontology")
    set_attribute(cv1, "URI",      "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo")
    set_attribute(cv1, "version",  "4.1.0")
    cv2 = new_child(cvList, "cv")
    set_attribute(cv2, "id",       "UO")
    set_attribute(cv2, "fullName", "Unit Ontology")
    set_attribute(cv2, "URI",      "http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo")
    set_attribute(cv2, "version",  "12:10:2011")
end

function _mzml_fileDescription(root::XMLElement)
    fd = new_child(root, "fileDescription")
    fc = new_child(fd, "fileContent")
    _cvParam(fc, "MS:1000579", "MS1 spectrum")
end

function _mzml_softwareList(root::XMLElement)
    sl = new_child(root, "softwareList")
    set_attribute(sl, "count", "1")
    sw = new_child(sl, "software")
    set_attribute(sw, "id",      "MassJ")
    set_attribute(sw, "version", "0.1")
    _cvParam(sw, "MS:1000799", "custom unreleased software tool"; value = "MassJ")
end

function _mzml_instrumentConfigurationList(root::XMLElement)
    icl = new_child(root, "instrumentConfigurationList")
    set_attribute(icl, "count", "1")
    ic = new_child(icl, "instrumentConfiguration")
    set_attribute(ic, "id", "IC1")
    _cvParam(ic, "MS:1000031", "instrument model")
end

function _mzml_dataProcessingList(root::XMLElement)
    dpl = new_child(root, "dataProcessingList")
    set_attribute(dpl, "count", "1")
    dp = new_child(dpl, "dataProcessing")
    set_attribute(dp, "id", "MassJExport")
    pm = new_child(dp, "processingMethod")
    set_attribute(pm, "order",       "0")
    set_attribute(pm, "softwareRef", "MassJ")
    _cvParam(pm, "MS:1000544", "Conversion to mzML")
end


# -- mzML spectrum body -------------------------------------------------------

function _mzml_spectrum(parent::XMLElement, scan::MSscan, index::Int;
                        precision::Int = 64, compress::Bool = true,
                        scalar::Bool = false)
    spec = new_child(parent, "spectrum")
    set_attribute(spec, "index",              string(index))
    set_attribute(spec, "id",                 "scan=$(scan.num)")
    set_attribute(spec, "defaultArrayLength", string(length(scan.mz)))

    if scalar
        sm = new_child(spec, "userParam")
        set_attribute(sm, "name",  MASSJ_SCALAR_PARAM)
        set_attribute(sm, "value", "true")
        set_attribute(sm, "type",  "xsd:string")
    end

    _cvParam(spec, CV_MS_LEVEL, "ms level"; value = string(scan.level))

    if scan.polarity == "+"
        _cvParam(spec, CV_POSITIVE_SCAN, "positive scan")
    elseif scan.polarity == "-"
        _cvParam(spec, CV_NEGATIVE_SCAN, "negative scan")
    end

    if scan.spectrumType === :centroid
        _cvParam(spec, CV_CENTROID, "centroid spectrum")
    elseif scan.spectrumType === :profile
        _cvParam(spec, CV_PROFILE, "profile spectrum")
    end

    _cvParam(spec, CV_TIC, "total ion current"; value = string(scan.tic))
    if scan.basePeakMz > 0
        _cvParam(spec, CV_BASE_PEAK_MZ, "base peak m/z";
                 value = string(scan.basePeakMz),
                 unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    end
    if scan.basePeakIntensity > 0
        _cvParam(spec, CV_BASE_PEAK_INT, "base peak intensity";
                 value = string(scan.basePeakIntensity),
                 unit_cv = "MS", unit_acc = "MS:1000131",
                 unit_name = "number of detector counts")
    end

    # Scan list with retention time
    scanList = new_child(spec, "scanList")
    set_attribute(scanList, "count", "1")
    _cvParam(scanList, "MS:1000795", "no combination")
    sc = new_child(scanList, "scan")
    _cvParam(sc, CV_SCAN_START_TIME, "scan start time";
             value     = string(scan.rt),
             unit_cv   = "UO",
             unit_acc  = CV_UNIT_MINUTE,
             unit_name = "minute")

    # Precursor list (for MS²+)
    if scan.level >= 2 && scan.precursor > 0
        precList = new_child(spec, "precursorList")
        set_attribute(precList, "count", "1")
        prec   = new_child(precList, "precursor")
        siList = new_child(prec, "selectedIonList")
        set_attribute(siList, "count", "1")
        si = new_child(siList, "selectedIon")
        _cvParam(si, CV_SELECTED_ION_MZ, "selected ion m/z";
                 value = string(scan.precursor),
                 unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        if scan.chargeState != 0
            _cvParam(si, CV_CHARGE_STATE, "charge state";
                     value = string(scan.chargeState))
        end
        act = new_child(prec, "activation")
        if !isempty(scan.activationMethod)
            for (accession, methodName) in ACTIVATION_METHODS
                if methodName == scan.activationMethod
                    _cvParam(act, accession, methodName)
                    break
                end
            end
        end
        if scan.collisionEnergy > 0
            _cvParam(act, CV_COLLISION_ENERGY, "collision energy";
                     value = string(scan.collisionEnergy),
                     unit_cv = "UO", unit_acc = "UO:0000266",
                     unit_name = "electronvolt")
        end
    end

    # Binary data arrays
    bdaList = new_child(spec, "binaryDataArrayList")
    set_attribute(bdaList, "count", "2")
    _mzml_binaryDataArray(bdaList, scan.mz,  :mz;  precision = precision, compress = compress)
    _mzml_binaryDataArray(bdaList, scan.int, :int; precision = precision, compress = compress)
end

function _mzml_binaryDataArray(parent::XMLElement, data::Vector{Float64},
                               kind::Symbol;
                               precision::Int = 64, compress::Bool = true)
    b64, _ = _encode_binary(data; precision = precision, compress = compress,
                            endian = :little)
    bda = new_child(parent, "binaryDataArray")
    set_attribute(bda, "encodedLength", string(length(b64)))

    prec_acc, prec_name = precision == 64 ?
        (CV_64BIT, "64-bit float") :
        (CV_32BIT, "32-bit float")
    _cvParam(bda, prec_acc, prec_name)

    comp_acc, comp_name = compress ?
        (CV_ZLIB, "zlib compression") :
        (CV_NO_COMPRESSION, "no compression")
    _cvParam(bda, comp_acc, comp_name)

    if kind === :mz
        _cvParam(bda, CV_MZ_ARRAY, "m/z array";
                 unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    else
        _cvParam(bda, CV_INT_ARRAY, "intensity array";
                 unit_cv = "MS", unit_acc = "MS:1000131",
                 unit_name = "number of detector counts")
    end

    bin = new_child(bda, "binary")
    add_text(bin, b64)
end


# Pipe-separator used to serialise vector-valued provenance fields into a
# single string-typed userParam (mzML) or custom attribute (mzXML).
const MASSJ_VEC_SEP = "|"

# Serialise an MSscans provenance vector as a `userParam` child of `parent`.
# Empty vectors are written as an empty `value`.
function _add_vec_userParam(parent::XMLElement, pname::String,
                            v::AbstractVector)
    up = new_child(parent, "userParam")
    set_attribute(up, "name",  pname)
    set_attribute(up, "value", isempty(v) ? "" : join(v, MASSJ_VEC_SEP))
    set_attribute(up, "type",  "xsd:string")
end

# Same idea for mzXML, which has no userParam — use a custom attribute.
function _set_vec_attr(elem::XMLElement, attrname::String, v::AbstractVector)
    set_attribute(elem, attrname, isempty(v) ? "" : join(v, MASSJ_VEC_SEP))
end


# -- mzML MSscans spectrum (with variance + marker) ---------------------------

function _mzml_msscans_spectrum(parent::XMLElement, scan::MSscans, index::Int = 0;
                                precision::Int = 64, compress::Bool = true,
                                scalar::Bool = true)
    spec = new_child(parent, "spectrum")
    num0 = isempty(scan.num) ? 1 : scan.num[1]
    set_attribute(spec, "index",              string(index))
    set_attribute(spec, "id",                 "scan=$(num0)")
    set_attribute(spec, "defaultArrayLength", string(length(scan.mz)))

    # Marker that this spectrum carries an averaged-spectrum payload
    up = new_child(spec, "userParam")
    set_attribute(up, "name",  MASSJ_CONTAINER_PARAM)
    set_attribute(up, "value", "MSscans")
    set_attribute(up, "type",  "xsd:string")

    # Scalar marker only when this MSscans was passed in bare, not as part of
    # a Vector{MSscans} (where load should keep the surrounding Vector).
    if scalar
        sm = new_child(spec, "userParam")
        set_attribute(sm, "name",  MASSJ_SCALAR_PARAM)
        set_attribute(sm, "value", "true")
        set_attribute(sm, "type",  "xsd:string")
    end

    # Preserve all vector-valued provenance fields so round-trip is loss-less.
    _add_vec_userParam(spec, "MassJ:num",                 scan.num)
    _add_vec_userParam(spec, "MassJ:rt",                  scan.rt)
    _add_vec_userParam(spec, "MassJ:level",               scan.level)
    _add_vec_userParam(spec, "MassJ:precursor",           scan.precursor)
    _add_vec_userParam(spec, "MassJ:polarity",            scan.polarity)
    _add_vec_userParam(spec, "MassJ:activationMethod",    scan.activationMethod)
    _add_vec_userParam(spec, "MassJ:collisionEnergy",     scan.collisionEnergy)
    _add_vec_userParam(spec, "MassJ:chargeState",         scan.chargeState)
    _add_vec_userParam(spec, "MassJ:driftTime",           scan.driftTime)
    _add_vec_userParam(spec, "MassJ:compensationVoltage", scan.compensationVoltage)

    lvl = isempty(scan.level) ? 1 : scan.level[1]
    _cvParam(spec, CV_MS_LEVEL, "ms level"; value = string(lvl))

    pol = isempty(scan.polarity) ? "" : scan.polarity[1]
    if pol == "+"
        _cvParam(spec, CV_POSITIVE_SCAN, "positive scan")
    elseif pol == "-"
        _cvParam(spec, CV_NEGATIVE_SCAN, "negative scan")
    end

    if scan.spectrumType === :centroid
        _cvParam(spec, CV_CENTROID, "centroid spectrum")
    elseif scan.spectrumType === :profile
        _cvParam(spec, CV_PROFILE, "profile spectrum")
    end

    _cvParam(spec, CV_TIC, "total ion current"; value = string(scan.tic))
    if scan.basePeakMz > 0
        _cvParam(spec, CV_BASE_PEAK_MZ, "base peak m/z";
                 value = string(scan.basePeakMz),
                 unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    end
    if scan.basePeakIntensity > 0
        _cvParam(spec, CV_BASE_PEAK_INT, "base peak intensity";
                 value = string(scan.basePeakIntensity),
                 unit_cv = "MS", unit_acc = "MS:1000131",
                 unit_name = "number of detector counts")
    end

    rt0 = isempty(scan.rt) ? 0.0 : scan.rt[1]
    scanList = new_child(spec, "scanList")
    set_attribute(scanList, "count", "1")
    _cvParam(scanList, "MS:1000795", "no combination")
    sc = new_child(scanList, "scan")
    _cvParam(sc, CV_SCAN_START_TIME, "scan start time";
             value = string(rt0), unit_cv = "UO",
             unit_acc = CV_UNIT_MINUTE, unit_name = "minute")

    prec0 = isempty(scan.precursor) ? 0.0 : scan.precursor[1]
    if lvl >= 2 && prec0 > 0
        precList = new_child(spec, "precursorList")
        set_attribute(precList, "count", "1")
        prec   = new_child(precList, "precursor")
        siList = new_child(prec, "selectedIonList")
        set_attribute(siList, "count", "1")
        si = new_child(siList, "selectedIon")
        _cvParam(si, CV_SELECTED_ION_MZ, "selected ion m/z";
                 value = string(prec0),
                 unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        chg0 = isempty(scan.chargeState) ? 0 : scan.chargeState[1]
        if chg0 != 0
            _cvParam(si, CV_CHARGE_STATE, "charge state"; value = string(chg0))
        end
        act = new_child(prec, "activation")
        am0 = isempty(scan.activationMethod) ? "" : scan.activationMethod[1]
        if !isempty(am0)
            for (accession, methodName) in ACTIVATION_METHODS
                if methodName == am0
                    _cvParam(act, accession, methodName)
                    break
                end
            end
        end
        ce0 = isempty(scan.collisionEnergy) ? 0.0 : scan.collisionEnergy[1]
        if ce0 > 0
            _cvParam(act, CV_COLLISION_ENERGY, "collision energy";
                     value = string(ce0),
                     unit_cv = "UO", unit_acc = "UO:0000266",
                     unit_name = "electronvolt")
        end
    end

    # Three binary arrays: m/z, intensity, variance (s)
    bdaList = new_child(spec, "binaryDataArrayList")
    set_attribute(bdaList, "count", "3")
    _mzml_binaryDataArray(bdaList, scan.mz,  :mz;  precision = precision, compress = compress)
    _mzml_binaryDataArray(bdaList, scan.int, :int; precision = precision, compress = compress)
    _mzml_binaryDataArray_variance(bdaList, scan.s;
                                   precision = precision, compress = compress)
end

function _mzml_binaryDataArray_variance(parent::XMLElement, data::Vector{Float64};
                                        precision::Int = 64, compress::Bool = true)
    b64, _ = _encode_binary(data; precision = precision, compress = compress,
                            endian = :little)
    bda = new_child(parent, "binaryDataArray")
    set_attribute(bda, "encodedLength", string(length(b64)))

    prec_acc, prec_name = precision == 64 ?
        (CV_64BIT, "64-bit float") : (CV_32BIT, "32-bit float")
    _cvParam(bda, prec_acc, prec_name)

    comp_acc, comp_name = compress ?
        (CV_ZLIB, "zlib compression") : (CV_NO_COMPRESSION, "no compression")
    _cvParam(bda, comp_acc, comp_name)

    # MassJ-specific marker — not a PSI-MS CV term, so emit as userParam.
    up = new_child(bda, "userParam")
    set_attribute(up, "name", MASSJ_VARIANCE_PARAM)
    set_attribute(up, "type", "xsd:string")

    bin = new_child(bda, "binary")
    add_text(bin, b64)
end


# ============================================================================
# mzXML writer
# ============================================================================

"""
    save_mzxml(filename::AbstractString, data;
               precision::Int = 64, compress::Bool = true) -> filename
Write a [`MSscan`](@ref), [`MSscans`](@ref), or `Vector{MSscan}` to an mzXML
file. The emitted file is minimal-but-valid and round-trips through
[`load`](@ref) — m/z, intensity, MS level, polarity, retention time, precursor
m/z, activation method, and collision energy are preserved. mzXML interleaves
m/z and intensity in a single `<peaks>` blob and uses *big-endian* byte order
("network"), in contrast with mzML.

Optional keywords:
* `precision = 64` — `64` for `Float64` arrays, `32` for `Float32`
* `compress = true` — zlib-compress the peaks blob
"""
function save_mzxml(filename::AbstractString, scans::Vector{MSscan};
                    precision::Int = 64, compress::Bool = true)
    return _save_mzxml_vector(filename, scans;
                              precision = precision, compress = compress,
                              scalar = false)
end

function save_mzxml(filename::AbstractString, scan::MSscan;
                    precision::Int = 64, compress::Bool = true)
    return _save_mzxml_vector(filename, [scan];
                              precision = precision, compress = compress,
                              scalar = true)
end

function _save_mzxml_vector(filename::AbstractString, scans::Vector{MSscan};
                            precision::Int, compress::Bool, scalar::Bool)
    xdoc  = XMLDocument()
    xroot = create_root(xdoc, "mzXML")
    set_attribute(xroot, "xmlns",     "http://sashimi.sourceforge.net/schema_revision/mzXML_3.2")
    set_attribute(xroot, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")

    msRun = new_child(xroot, "msRun")
    set_attribute(msRun, "scanCount", string(length(scans)))

    for scan in scans
        _mzxml_spectrum(msRun, scan;
                        precision = precision, compress = compress,
                        scalar = scalar)
    end

    save_file(xdoc, filename)
    free(xdoc)
    return filename
end

function save_mzxml(filename::AbstractString, scan::MSscans;
                    precision::Int = 64, compress::Bool = true)
    return _save_mzxml_msscans_vector(filename, [scan];
                                      precision = precision, compress = compress,
                                      scalar = true)
end

function save_mzxml(filename::AbstractString, scans::Vector{MSscans};
                    precision::Int = 64, compress::Bool = true)
    return _save_mzxml_msscans_vector(filename, scans;
                                      precision = precision, compress = compress,
                                      scalar = false)
end

function _save_mzxml_msscans_vector(filename::AbstractString, scans::Vector{MSscans};
                                    precision::Int, compress::Bool, scalar::Bool)
    xdoc  = XMLDocument()
    xroot = create_root(xdoc, "mzXML")
    set_attribute(xroot, "xmlns",     "http://sashimi.sourceforge.net/schema_revision/mzXML_3.2")
    set_attribute(xroot, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")

    msRun = new_child(xroot, "msRun")
    set_attribute(msRun, "scanCount", string(length(scans)))

    for sc in scans
        _mzxml_msscans_spectrum(msRun, sc;
                                precision = precision, compress = compress,
                                scalar = scalar)
    end

    save_file(xdoc, filename)
    free(xdoc)
    return filename
end


function _mzxml_msscans_spectrum(parent::XMLElement, scan::MSscans;
                                 precision::Int = 64, compress::Bool = true,
                                 scalar::Bool = true)
    sc = new_child(parent, "scan")
    num0 = isempty(scan.num) ? 1 : scan.num[1]
    set_attribute(sc, "num",        string(num0))
    lvl  = isempty(scan.level) ? 1 : scan.level[1]
    set_attribute(sc, "msLevel",    string(lvl))
    set_attribute(sc, "peaksCount", string(length(scan.mz)))
    # MassJ-specific markers: averaged container + optional saved-as-scalar
    set_attribute(sc, MASSJ_MZXML_CONTAINER_ATTR, "MSscans")
    if scalar
        set_attribute(sc, MASSJ_MZXML_SCALAR_ATTR, "true")
    end

    # Vector-valued provenance fields, joined with | as a custom attribute each.
    _set_vec_attr(sc, "MassJNum",                 scan.num)
    _set_vec_attr(sc, "MassJRt",                  scan.rt)
    _set_vec_attr(sc, "MassJLevel",               scan.level)
    _set_vec_attr(sc, "MassJPrecursor",           scan.precursor)
    _set_vec_attr(sc, "MassJPolarity",            scan.polarity)
    _set_vec_attr(sc, "MassJActivationMethod",    scan.activationMethod)
    _set_vec_attr(sc, "MassJCollisionEnergy",     scan.collisionEnergy)
    _set_vec_attr(sc, "MassJChargeState",         scan.chargeState)
    _set_vec_attr(sc, "MassJDriftTime",           scan.driftTime)
    _set_vec_attr(sc, "MassJCompensationVoltage", scan.compensationVoltage)

    pol = isempty(scan.polarity) ? "" : scan.polarity[1]
    if !isempty(pol)
        set_attribute(sc, "polarity", pol)
    end
    rt0 = isempty(scan.rt) ? 0.0 : scan.rt[1]
    set_attribute(sc, "retentionTime", "PT$(rt0)M")
    set_attribute(sc, "totIonCurrent", string(scan.tic))
    if scan.basePeakMz > 0
        set_attribute(sc, "basePeakMz",        string(scan.basePeakMz))
        set_attribute(sc, "basePeakIntensity", string(scan.basePeakIntensity))
    end
    ce0 = isempty(scan.collisionEnergy) ? 0.0 : scan.collisionEnergy[1]
    if ce0 > 0
        set_attribute(sc, "collisionEnergy", string(ce0))
    end

    prec0 = isempty(scan.precursor) ? 0.0 : scan.precursor[1]
    if lvl >= 2 && prec0 > 0
        pm = new_child(sc, "precursorMz")
        am0 = isempty(scan.activationMethod) ? "" : scan.activationMethod[1]
        if !isempty(am0)
            set_attribute(pm, "activationMethod", am0)
        end
        add_text(pm, string(prec0))
    end

    # Standard interleaved (m/z, intensity) peaks blob
    n = length(scan.mz)
    interleaved = Vector{Float64}(undef, 2n)
    @inbounds for i in 1:n
        interleaved[2i - 1] = scan.mz[i]
        interleaved[2i]     = scan.int[i]
    end
    b64, byte_len = _encode_binary(interleaved;
                                   precision = precision,
                                   compress  = compress,
                                   endian    = :big)
    peaks = new_child(sc, "peaks")
    set_attribute(peaks, "precision",       string(precision))
    set_attribute(peaks, "byteOrder",       "network")
    set_attribute(peaks, "pairOrder",       "m/z-int")
    set_attribute(peaks, "contentType",     "m/z-int")
    set_attribute(peaks, "compressionType", compress ? "zlib" : "none")
    compress && set_attribute(peaks, "compressedLen", string(byte_len))
    add_text(peaks, b64)

    # Variance blob — a second <peaks> child with pairOrder="variance"
    b64v, byte_len_v = _encode_binary(scan.s;
                                      precision = precision,
                                      compress  = compress,
                                      endian    = :big)
    vpeaks = new_child(sc, "peaks")
    set_attribute(vpeaks, "precision",       string(precision))
    set_attribute(vpeaks, "byteOrder",       "network")
    set_attribute(vpeaks, "pairOrder",       MASSJ_MZXML_VARIANCE_PAIR)
    set_attribute(vpeaks, "contentType",     MASSJ_MZXML_VARIANCE_PAIR)
    set_attribute(vpeaks, "compressionType", compress ? "zlib" : "none")
    compress && set_attribute(vpeaks, "compressedLen", string(byte_len_v))
    add_text(vpeaks, b64v)
end


function _mzxml_spectrum(parent::XMLElement, scan::MSscan;
                         precision::Int = 64, compress::Bool = true,
                         scalar::Bool = false)
    sc = new_child(parent, "scan")
    set_attribute(sc, "num",     string(scan.num))
    set_attribute(sc, "msLevel", string(scan.level))
    set_attribute(sc, "peaksCount", string(length(scan.mz)))
    if scalar
        set_attribute(sc, MASSJ_MZXML_SCALAR_ATTR, "true")
    end
    if !isempty(scan.polarity)
        set_attribute(sc, "polarity", scan.polarity)
    end
    # ISO-8601 duration in minutes, e.g. "PT0.1384M" — read side strips "PT…M"
    set_attribute(sc, "retentionTime", "PT$(scan.rt)M")
    set_attribute(sc, "totIonCurrent", string(scan.tic))
    if scan.basePeakMz > 0
        set_attribute(sc, "basePeakMz",        string(scan.basePeakMz))
        set_attribute(sc, "basePeakIntensity", string(scan.basePeakIntensity))
    end
    if scan.collisionEnergy > 0
        set_attribute(sc, "collisionEnergy", string(scan.collisionEnergy))
    end

    if scan.level >= 2 && scan.precursor > 0
        pm = new_child(sc, "precursorMz")
        if !isempty(scan.activationMethod)
            set_attribute(pm, "activationMethod", scan.activationMethod)
        end
        add_text(pm, string(scan.precursor))
    end

    # Interleave m/z and intensity: [mz1, int1, mz2, int2, ...]
    n = length(scan.mz)
    interleaved = Vector{Float64}(undef, 2n)
    @inbounds for i in 1:n
        interleaved[2i - 1] = scan.mz[i]
        interleaved[2i]     = scan.int[i]
    end

    b64, byte_len = _encode_binary(interleaved;
                                   precision = precision,
                                   compress  = compress,
                                   endian    = :big)

    peaks = new_child(sc, "peaks")
    set_attribute(peaks, "precision",     string(precision))
    set_attribute(peaks, "byteOrder",     "network")
    set_attribute(peaks, "pairOrder",     "m/z-int")
    set_attribute(peaks, "contentType",   "m/z-int")
    set_attribute(peaks, "compressionType", compress ? "zlib" : "none")
    if compress
        set_attribute(peaks, "compressedLen", string(byte_len))
    end
    add_text(peaks, b64)
end
