"""
Export MassJ data to standard mass-spectrometry file formats.

The top-level entry point [`save`](@ref) dispatches on the file extension:

* `.mzML`  → [`save_mzml`](@ref)
* `.mzXML` → [`save_mzxml`](@ref)

Both writers stream directly to disk, one spectrum at a time, so peak memory
stays bounded by the largest single spectrum (not by the total file size). The
output round-trips through MassJ's own readers — `load(save_path)` recovers
the same data. Format-level metadata (`fileDescription`, `instrumentConfiguration`,
`dataProcessing`, the `cvList`) is emitted in a minimal-but-valid form.
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

Accepts a single [`MSscan`](@ref) / [`MSscans`](@ref), a `Vector{MSscan}`, or
a `Vector{MSscans}`. Common keyword arguments:

* `precision::Int = 64`   — `64` for `Float64` arrays, `32` for `Float32`
* `compress::Bool = true` — zlib-compress the binary arrays
* `progress::Bool = true` — show a [`ProgressMeter`](https://github.com/timholy/ProgressMeter.jl)
  progress bar while writing (set `false` in scripts / CI)

# Examples
```julia
scans = load("input.mzML")
save(scans, "output.mzML")                          # round-trip
save(scans, "output.mzXML"; precision = 32)         # smaller file, lossy
save(scans[1], "single_scan.mzML")                  # one scan
save(scans, "quiet.mzML"; progress = false)         # no progress bar
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


# ----------------------------------------------------------------------------
# Binary encoding (shared)
# ----------------------------------------------------------------------------

"""
    _encode_binary(data; precision = 64, compress = true, endian = :little)
        -> (base64_string, byte_length)
Encode a real-valued vector to the base64-of-zlib-of-raw-bytes payload expected
by mzML/mzXML `<binary>` / `<peaks>` elements. `endian` is `:little` for mzML
or `:big` for mzXML (network byte order).
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


# ----------------------------------------------------------------------------
# XML helpers — write directly to an IO, no in-memory tree.
# ----------------------------------------------------------------------------

# Escape the five XML metacharacters in attribute values / text content. Most
# values we emit are numeric and don't need escaping; user-supplied labels
# might.
function _xmlescape(s::AbstractString)
    needs = false
    for c in s
        if c == '&' || c == '<' || c == '>' || c == '"' || c == '\''
            needs = true
            break
        end
    end
    needs || return s
    buf = IOBuffer()
    for c in s
        c == '&'  ? print(buf, "&amp;")  :
        c == '<'  ? print(buf, "&lt;")   :
        c == '>'  ? print(buf, "&gt;")   :
        c == '"'  ? print(buf, "&quot;") :
        c == '\'' ? print(buf, "&apos;") :
                    print(buf, c)
    end
    return String(take!(buf))
end


# Emit `<cvParam cvRef="MS" accession="..." name="..." [value="..."] [unit*]/>`
function _stream_cvParam(io::IO, accession::AbstractString, pname::AbstractString;
                        value::AbstractString    = "",
                        unit_cv::AbstractString  = "",
                        unit_acc::AbstractString = "",
                        unit_name::AbstractString = "")
    write(io, "<cvParam cvRef=\"MS\" accession=\"", accession,
          "\" name=\"", _xmlescape(pname), "\"")
    if !isempty(value)
        write(io, " value=\"", _xmlescape(value), "\"")
    end
    if !isempty(unit_acc)
        write(io, " unitCvRef=\"", unit_cv,
              "\" unitAccession=\"", unit_acc,
              "\" unitName=\"", _xmlescape(unit_name), "\"")
    end
    write(io, "/>\n")
end


# Emit `<userParam name="..." [value="..."] type="xsd:string"/>`
function _stream_userParam(io::IO, pname::AbstractString;
                          value::AbstractString = "")
    write(io, "<userParam name=\"", _xmlescape(pname), "\"")
    if !isempty(value)
        write(io, " value=\"", _xmlescape(value), "\"")
    end
    write(io, " type=\"xsd:string\"/>\n")
end


# Serialise an MSscans provenance vector as a `userParam` (mzML).
function _stream_vec_userParam(io::IO, pname::AbstractString, v::AbstractVector)
    _stream_userParam(io, pname; value = isempty(v) ? "" : join(v, "|"))
end


# ----------------------------------------------------------------------------
# Streaming mzML writer
# ----------------------------------------------------------------------------

function _stream_mzml_open(io::IO, scancount::Int)
    write(io, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
    write(io, "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"",
          " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
          " version=\"1.1.0\">\n")

    write(io, "<cvList count=\"2\">\n",
              "<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\"",
              " URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"",
              " version=\"4.1.0\"/>\n",
              "<cv id=\"UO\" fullName=\"Unit Ontology\"",
              " URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"",
              " version=\"12:10:2011\"/>\n",
              "</cvList>\n")

    write(io, "<fileDescription>\n<fileContent>\n")
    _stream_cvParam(io, "MS:1000579", "MS1 spectrum")
    write(io, "</fileContent>\n</fileDescription>\n")

    write(io, "<softwareList count=\"1\">\n",
              "<software id=\"MassJ\" version=\"0.1\">\n")
    _stream_cvParam(io, "MS:1000799", "custom unreleased software tool"; value = "MassJ")
    write(io, "</software>\n</softwareList>\n")

    write(io, "<instrumentConfigurationList count=\"1\">\n",
              "<instrumentConfiguration id=\"IC1\">\n")
    _stream_cvParam(io, "MS:1000031", "instrument model")
    write(io, "</instrumentConfiguration>\n</instrumentConfigurationList>\n")

    write(io, "<dataProcessingList count=\"1\">\n",
              "<dataProcessing id=\"MassJExport\">\n",
              "<processingMethod order=\"0\" softwareRef=\"MassJ\">\n")
    _stream_cvParam(io, "MS:1000544", "Conversion to mzML")
    write(io, "</processingMethod>\n</dataProcessing>\n</dataProcessingList>\n")

    write(io, "<run id=\"run1\" defaultInstrumentConfigurationRef=\"IC1\">\n",
              "<spectrumList count=\"", string(scancount),
              "\" defaultDataProcessingRef=\"MassJExport\">\n")
end


_stream_mzml_close(io::IO) = write(io, "</spectrumList>\n</run>\n</mzML>\n")


function _stream_mzml_spectrum(io::IO, scan::MSscan, index::Int;
                               precision::Int = 64, compress::Bool = true,
                               scalar::Bool = false)
    write(io, "<spectrum index=\"", string(index),
              "\" id=\"scan=", string(scan.num),
              "\" defaultArrayLength=\"", string(length(scan.mz)), "\">\n")

    if scalar
        _stream_userParam(io, MASSJ_SCALAR_PARAM; value = "true")
    end

    _stream_cvParam(io, CV_MS_LEVEL, "ms level"; value = string(scan.level))

    if scan.polarity == "+"
        _stream_cvParam(io, CV_POSITIVE_SCAN, "positive scan")
    elseif scan.polarity == "-"
        _stream_cvParam(io, CV_NEGATIVE_SCAN, "negative scan")
    end

    if scan.spectrumType === :centroid
        _stream_cvParam(io, CV_CENTROID, "centroid spectrum")
    elseif scan.spectrumType === :profile
        _stream_cvParam(io, CV_PROFILE, "profile spectrum")
    end

    _stream_cvParam(io, CV_TIC, "total ion current"; value = string(scan.tic))
    if scan.basePeakMz > 0
        _stream_cvParam(io, CV_BASE_PEAK_MZ, "base peak m/z";
                       value = string(scan.basePeakMz),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    end
    if scan.basePeakIntensity > 0
        _stream_cvParam(io, CV_BASE_PEAK_INT, "base peak intensity";
                       value = string(scan.basePeakIntensity),
                       unit_cv = "MS", unit_acc = "MS:1000131",
                       unit_name = "number of detector counts")
    end

    write(io, "<scanList count=\"1\">\n")
    _stream_cvParam(io, "MS:1000795", "no combination")
    write(io, "<scan>\n")
    _stream_cvParam(io, CV_SCAN_START_TIME, "scan start time";
                   value = string(scan.rt), unit_cv = "UO",
                   unit_acc = CV_UNIT_MINUTE, unit_name = "minute")
    write(io, "</scan>\n</scanList>\n")

    if scan.level >= 2 && scan.precursor > 0
        write(io, "<precursorList count=\"1\">\n<precursor>\n",
                  "<selectedIonList count=\"1\">\n<selectedIon>\n")
        _stream_cvParam(io, CV_SELECTED_ION_MZ, "selected ion m/z";
                       value = string(scan.precursor),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        if scan.chargeState != 0
            _stream_cvParam(io, CV_CHARGE_STATE, "charge state";
                           value = string(scan.chargeState))
        end
        write(io, "</selectedIon>\n</selectedIonList>\n<activation>\n")
        if !isempty(scan.activationMethod)
            for (accession, methodName) in ACTIVATION_METHODS
                if methodName == scan.activationMethod
                    _stream_cvParam(io, accession, methodName)
                    break
                end
            end
        end
        if scan.collisionEnergy > 0
            _stream_cvParam(io, CV_COLLISION_ENERGY, "collision energy";
                           value = string(scan.collisionEnergy),
                           unit_cv = "UO", unit_acc = "UO:0000266",
                           unit_name = "electronvolt")
        end
        write(io, "</activation>\n</precursor>\n</precursorList>\n")
    end

    write(io, "<binaryDataArrayList count=\"2\">\n")
    _stream_mzml_binaryDataArray(io, scan.mz,  :mz;  precision = precision, compress = compress)
    _stream_mzml_binaryDataArray(io, scan.int, :int; precision = precision, compress = compress)
    write(io, "</binaryDataArrayList>\n</spectrum>\n")
end


function _stream_mzml_msscans_spectrum(io::IO, scan::MSscans, index::Int;
                                       precision::Int = 64, compress::Bool = true,
                                       scalar::Bool = true)
    num0 = isempty(scan.num) ? 1 : scan.num[1]
    write(io, "<spectrum index=\"", string(index),
              "\" id=\"scan=", string(num0),
              "\" defaultArrayLength=\"", string(length(scan.mz)), "\">\n")

    _stream_userParam(io, MASSJ_CONTAINER_PARAM; value = "MSscans")
    if scalar
        _stream_userParam(io, MASSJ_SCALAR_PARAM; value = "true")
    end

    # Preserve all vector-valued provenance fields losslessly.
    _stream_vec_userParam(io, "MassJ:num",                 scan.num)
    _stream_vec_userParam(io, "MassJ:rt",                  scan.rt)
    _stream_vec_userParam(io, "MassJ:level",               scan.level)
    _stream_vec_userParam(io, "MassJ:precursor",           scan.precursor)
    _stream_vec_userParam(io, "MassJ:polarity",            scan.polarity)
    _stream_vec_userParam(io, "MassJ:activationMethod",    scan.activationMethod)
    _stream_vec_userParam(io, "MassJ:collisionEnergy",     scan.collisionEnergy)
    _stream_vec_userParam(io, "MassJ:chargeState",         scan.chargeState)
    _stream_vec_userParam(io, "MassJ:driftTime",           scan.driftTime)
    _stream_vec_userParam(io, "MassJ:compensationVoltage", scan.compensationVoltage)

    lvl = isempty(scan.level) ? 1 : scan.level[1]
    _stream_cvParam(io, CV_MS_LEVEL, "ms level"; value = string(lvl))

    pol = isempty(scan.polarity) ? "" : scan.polarity[1]
    if pol == "+"
        _stream_cvParam(io, CV_POSITIVE_SCAN, "positive scan")
    elseif pol == "-"
        _stream_cvParam(io, CV_NEGATIVE_SCAN, "negative scan")
    end

    if scan.spectrumType === :centroid
        _stream_cvParam(io, CV_CENTROID, "centroid spectrum")
    elseif scan.spectrumType === :profile
        _stream_cvParam(io, CV_PROFILE, "profile spectrum")
    end

    _stream_cvParam(io, CV_TIC, "total ion current"; value = string(scan.tic))
    if scan.basePeakMz > 0
        _stream_cvParam(io, CV_BASE_PEAK_MZ, "base peak m/z";
                       value = string(scan.basePeakMz),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    end
    if scan.basePeakIntensity > 0
        _stream_cvParam(io, CV_BASE_PEAK_INT, "base peak intensity";
                       value = string(scan.basePeakIntensity),
                       unit_cv = "MS", unit_acc = "MS:1000131",
                       unit_name = "number of detector counts")
    end

    rt0 = isempty(scan.rt) ? 0.0 : scan.rt[1]
    write(io, "<scanList count=\"1\">\n")
    _stream_cvParam(io, "MS:1000795", "no combination")
    write(io, "<scan>\n")
    _stream_cvParam(io, CV_SCAN_START_TIME, "scan start time";
                   value = string(rt0), unit_cv = "UO",
                   unit_acc = CV_UNIT_MINUTE, unit_name = "minute")
    write(io, "</scan>\n</scanList>\n")

    prec0 = isempty(scan.precursor) ? 0.0 : scan.precursor[1]
    if lvl >= 2 && prec0 > 0
        write(io, "<precursorList count=\"1\">\n<precursor>\n",
                  "<selectedIonList count=\"1\">\n<selectedIon>\n")
        _stream_cvParam(io, CV_SELECTED_ION_MZ, "selected ion m/z";
                       value = string(prec0),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        chg0 = isempty(scan.chargeState) ? 0 : scan.chargeState[1]
        if chg0 != 0
            _stream_cvParam(io, CV_CHARGE_STATE, "charge state"; value = string(chg0))
        end
        write(io, "</selectedIon>\n</selectedIonList>\n<activation>\n")
        am0 = isempty(scan.activationMethod) ? "" : scan.activationMethod[1]
        if !isempty(am0)
            for (accession, methodName) in ACTIVATION_METHODS
                if methodName == am0
                    _stream_cvParam(io, accession, methodName)
                    break
                end
            end
        end
        ce0 = isempty(scan.collisionEnergy) ? 0.0 : scan.collisionEnergy[1]
        if ce0 > 0
            _stream_cvParam(io, CV_COLLISION_ENERGY, "collision energy";
                           value = string(ce0),
                           unit_cv = "UO", unit_acc = "UO:0000266",
                           unit_name = "electronvolt")
        end
        write(io, "</activation>\n</precursor>\n</precursorList>\n")
    end

    write(io, "<binaryDataArrayList count=\"3\">\n")
    _stream_mzml_binaryDataArray(io, scan.mz,  :mz;  precision = precision, compress = compress)
    _stream_mzml_binaryDataArray(io, scan.int, :int; precision = precision, compress = compress)
    _stream_mzml_variance_array(io, scan.s; precision = precision, compress = compress)
    write(io, "</binaryDataArrayList>\n</spectrum>\n")
end


function _stream_mzml_binaryDataArray(io::IO, data::Vector{Float64}, kind::Symbol;
                                      precision::Int = 64, compress::Bool = true)
    b64, _ = _encode_binary(data; precision = precision, compress = compress,
                            endian = :little)
    write(io, "<binaryDataArray encodedLength=\"", string(length(b64)), "\">\n")

    if precision == 64
        _stream_cvParam(io, CV_64BIT, "64-bit float")
    else
        _stream_cvParam(io, CV_32BIT, "32-bit float")
    end
    if compress
        _stream_cvParam(io, CV_ZLIB, "zlib compression")
    else
        _stream_cvParam(io, CV_NO_COMPRESSION, "no compression")
    end
    if kind === :mz
        _stream_cvParam(io, CV_MZ_ARRAY, "m/z array";
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
    else
        _stream_cvParam(io, CV_INT_ARRAY, "intensity array";
                       unit_cv = "MS", unit_acc = "MS:1000131",
                       unit_name = "number of detector counts")
    end

    write(io, "<binary>", b64, "</binary>\n</binaryDataArray>\n")
end


function _stream_mzml_variance_array(io::IO, data::Vector{Float64};
                                     precision::Int = 64, compress::Bool = true)
    b64, _ = _encode_binary(data; precision = precision, compress = compress,
                            endian = :little)
    write(io, "<binaryDataArray encodedLength=\"", string(length(b64)), "\">\n")

    if precision == 64
        _stream_cvParam(io, CV_64BIT, "64-bit float")
    else
        _stream_cvParam(io, CV_32BIT, "32-bit float")
    end
    if compress
        _stream_cvParam(io, CV_ZLIB, "zlib compression")
    else
        _stream_cvParam(io, CV_NO_COMPRESSION, "no compression")
    end
    _stream_userParam(io, MASSJ_VARIANCE_PARAM)
    write(io, "<binary>", b64, "</binary>\n</binaryDataArray>\n")
end


# -- save_mzml dispatch ------------------------------------------------------

"""
    save_mzml(filename::AbstractString, data;
              precision::Int = 64, compress::Bool = true,
              progress::Bool = true) -> filename
Write a [`MSscan`](@ref), [`MSscans`](@ref), or `Vector{MSscan}` /
`Vector{MSscans}` to an mzML file. The file is written one spectrum at a time;
peak RAM is bounded by the largest single spectrum, not the total file size.

When `progress = true` (default) a `ProgressMeter` bar is shown while writing.

Round-trips through [`load`](@ref): the loaded value has the same type as the
saved one (scalar `MSscan` / `MSscans` come back bare; vectors come back as
vectors).
"""
function save_mzml(filename::AbstractString, scans::Vector{MSscan};
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true)
    return _save_mzml_vector(filename, scans;
                             precision = precision, compress = compress,
                             scalar = false, progress = progress)
end

function save_mzml(filename::AbstractString, scan::MSscan;
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true)
    return _save_mzml_vector(filename, [scan];
                             precision = precision, compress = compress,
                             scalar = true, progress = progress)
end

function _save_mzml_vector(filename::AbstractString, scans::Vector{MSscan};
                           precision::Int, compress::Bool,
                           scalar::Bool, progress::Bool)
    open(filename, "w") do io
        _stream_mzml_open(io, length(scans))
        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzML: ") : nothing
        for (i, scan) in enumerate(scans)
            _stream_mzml_spectrum(io, scan, i - 1;
                                  precision = precision, compress = compress,
                                  scalar = scalar)
            prog === nothing || next!(prog)
        end
        _stream_mzml_close(io)
    end
    return filename
end


function save_mzml(filename::AbstractString, scan::MSscans;
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true)
    return _save_mzml_msscans_vector(filename, [scan];
                                     precision = precision, compress = compress,
                                     scalar = true, progress = progress)
end

function save_mzml(filename::AbstractString, scans::Vector{MSscans};
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true)
    return _save_mzml_msscans_vector(filename, scans;
                                     precision = precision, compress = compress,
                                     scalar = false, progress = progress)
end

function _save_mzml_msscans_vector(filename::AbstractString, scans::Vector{MSscans};
                                   precision::Int, compress::Bool,
                                   scalar::Bool, progress::Bool)
    open(filename, "w") do io
        _stream_mzml_open(io, length(scans))
        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzML: ") : nothing
        for (i, sc) in enumerate(scans)
            _stream_mzml_msscans_spectrum(io, sc, i - 1;
                                          precision = precision, compress = compress,
                                          scalar = scalar)
            prog === nothing || next!(prog)
        end
        _stream_mzml_close(io)
    end
    return filename
end


# ----------------------------------------------------------------------------
# Streaming mzXML writer
# ----------------------------------------------------------------------------

function _stream_mzxml_open(io::IO, scancount::Int)
    write(io, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
              "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2\"",
              " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
              "<msRun scanCount=\"", string(scancount), "\">\n")
end


_stream_mzxml_close(io::IO) = write(io, "</msRun>\n</mzXML>\n")


# Pipe-joined vector attribute for mzXML.
function _attr_vec(v::AbstractVector)
    isempty(v) ? "" : join(v, "|")
end


function _stream_mzxml_spectrum(io::IO, scan::MSscan;
                                precision::Int = 64, compress::Bool = true,
                                scalar::Bool = false)
    write(io, "<scan num=\"", string(scan.num),
              "\" msLevel=\"", string(scan.level),
              "\" peaksCount=\"", string(length(scan.mz)), "\"")
    if scalar
        write(io, " ", MASSJ_MZXML_SCALAR_ATTR, "=\"true\"")
    end
    if !isempty(scan.polarity)
        write(io, " polarity=\"", _xmlescape(scan.polarity), "\"")
    end
    write(io, " retentionTime=\"PT", string(scan.rt), "M\"")
    write(io, " totIonCurrent=\"", string(scan.tic), "\"")
    if scan.basePeakMz > 0
        write(io, " basePeakMz=\"",        string(scan.basePeakMz),
                  "\" basePeakIntensity=\"", string(scan.basePeakIntensity), "\"")
    end
    if scan.collisionEnergy > 0
        write(io, " collisionEnergy=\"", string(scan.collisionEnergy), "\"")
    end
    write(io, ">\n")

    if scan.level >= 2 && scan.precursor > 0
        write(io, "<precursorMz")
        if !isempty(scan.activationMethod)
            write(io, " activationMethod=\"", _xmlescape(scan.activationMethod), "\"")
        end
        write(io, ">", string(scan.precursor), "</precursorMz>\n")
    end

    # Interleaved m/z, intensity, m/z, intensity, …
    n = length(scan.mz)
    interleaved = Vector{Float64}(undef, 2n)
    @inbounds for i in 1:n
        interleaved[2i - 1] = scan.mz[i]
        interleaved[2i]     = scan.int[i]
    end
    b64, byte_len = _encode_binary(interleaved;
                                   precision = precision, compress = compress,
                                   endian = :big)
    write(io, "<peaks precision=\"", string(precision),
              "\" byteOrder=\"network\" pairOrder=\"m/z-int\" contentType=\"m/z-int\"",
              " compressionType=\"", compress ? "zlib" : "none", "\"")
    if compress
        write(io, " compressedLen=\"", string(byte_len), "\"")
    end
    write(io, ">", b64, "</peaks>\n")

    write(io, "</scan>\n")
end


function _stream_mzxml_msscans_spectrum(io::IO, scan::MSscans;
                                        precision::Int = 64, compress::Bool = true,
                                        scalar::Bool = true)
    num0 = isempty(scan.num) ? 1 : scan.num[1]
    lvl  = isempty(scan.level) ? 1 : scan.level[1]

    write(io, "<scan num=\"", string(num0),
              "\" msLevel=\"", string(lvl),
              "\" peaksCount=\"", string(length(scan.mz)), "\"")

    # MassJ markers + serialised vector provenance, all as custom attributes.
    write(io, " ", MASSJ_MZXML_CONTAINER_ATTR, "=\"MSscans\"")
    if scalar
        write(io, " ", MASSJ_MZXML_SCALAR_ATTR, "=\"true\"")
    end
    write(io, " MassJNum=\"",                 _xmlescape(_attr_vec(scan.num)),                 "\"")
    write(io, " MassJRt=\"",                  _xmlescape(_attr_vec(scan.rt)),                  "\"")
    write(io, " MassJLevel=\"",               _xmlescape(_attr_vec(scan.level)),               "\"")
    write(io, " MassJPrecursor=\"",           _xmlescape(_attr_vec(scan.precursor)),           "\"")
    write(io, " MassJPolarity=\"",            _xmlescape(_attr_vec(scan.polarity)),            "\"")
    write(io, " MassJActivationMethod=\"",    _xmlescape(_attr_vec(scan.activationMethod)),    "\"")
    write(io, " MassJCollisionEnergy=\"",     _xmlescape(_attr_vec(scan.collisionEnergy)),     "\"")
    write(io, " MassJChargeState=\"",         _xmlescape(_attr_vec(scan.chargeState)),         "\"")
    write(io, " MassJDriftTime=\"",           _xmlescape(_attr_vec(scan.driftTime)),           "\"")
    write(io, " MassJCompensationVoltage=\"", _xmlescape(_attr_vec(scan.compensationVoltage)), "\"")

    pol = isempty(scan.polarity) ? "" : scan.polarity[1]
    if !isempty(pol)
        write(io, " polarity=\"", _xmlescape(pol), "\"")
    end
    rt0 = isempty(scan.rt) ? 0.0 : scan.rt[1]
    write(io, " retentionTime=\"PT", string(rt0), "M\"")
    write(io, " totIonCurrent=\"", string(scan.tic), "\"")
    if scan.basePeakMz > 0
        write(io, " basePeakMz=\"",        string(scan.basePeakMz),
                  "\" basePeakIntensity=\"", string(scan.basePeakIntensity), "\"")
    end
    ce0 = isempty(scan.collisionEnergy) ? 0.0 : scan.collisionEnergy[1]
    if ce0 > 0
        write(io, " collisionEnergy=\"", string(ce0), "\"")
    end
    write(io, ">\n")

    prec0 = isempty(scan.precursor) ? 0.0 : scan.precursor[1]
    if lvl >= 2 && prec0 > 0
        write(io, "<precursorMz")
        am0 = isempty(scan.activationMethod) ? "" : scan.activationMethod[1]
        if !isempty(am0)
            write(io, " activationMethod=\"", _xmlescape(am0), "\"")
        end
        write(io, ">", string(prec0), "</precursorMz>\n")
    end

    # Standard m/z+intensity peaks blob
    n = length(scan.mz)
    interleaved = Vector{Float64}(undef, 2n)
    @inbounds for i in 1:n
        interleaved[2i - 1] = scan.mz[i]
        interleaved[2i]     = scan.int[i]
    end
    b64, byte_len = _encode_binary(interleaved;
                                   precision = precision, compress = compress,
                                   endian = :big)
    write(io, "<peaks precision=\"", string(precision),
              "\" byteOrder=\"network\" pairOrder=\"m/z-int\" contentType=\"m/z-int\"",
              " compressionType=\"", compress ? "zlib" : "none", "\"")
    if compress
        write(io, " compressedLen=\"", string(byte_len), "\"")
    end
    write(io, ">", b64, "</peaks>\n")

    # Variance blob — second <peaks> child
    b64v, byte_len_v = _encode_binary(scan.s;
                                       precision = precision, compress = compress,
                                       endian = :big)
    write(io, "<peaks precision=\"", string(precision),
              "\" byteOrder=\"network\" pairOrder=\"", MASSJ_MZXML_VARIANCE_PAIR,
              "\" contentType=\"", MASSJ_MZXML_VARIANCE_PAIR, "\"",
              " compressionType=\"", compress ? "zlib" : "none", "\"")
    if compress
        write(io, " compressedLen=\"", string(byte_len_v), "\"")
    end
    write(io, ">", b64v, "</peaks>\n")

    write(io, "</scan>\n")
end


# -- save_mzxml dispatch -----------------------------------------------------

"""
    save_mzxml(filename::AbstractString, data;
               precision::Int = 64, compress::Bool = true,
               progress::Bool = true) -> filename
Write a [`MSscan`](@ref), [`MSscans`](@ref), or `Vector{MSscan}` /
`Vector{MSscans}` to an mzXML file. Streams one spectrum at a time; peak RAM is
bounded by the largest single spectrum. mzXML interleaves m/z and intensity in
a single `<peaks>` blob and uses big-endian (network) byte order.
"""
function save_mzxml(filename::AbstractString, scans::Vector{MSscan};
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true)
    return _save_mzxml_vector(filename, scans;
                              precision = precision, compress = compress,
                              scalar = false, progress = progress)
end

function save_mzxml(filename::AbstractString, scan::MSscan;
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true)
    return _save_mzxml_vector(filename, [scan];
                              precision = precision, compress = compress,
                              scalar = true, progress = progress)
end

function _save_mzxml_vector(filename::AbstractString, scans::Vector{MSscan};
                            precision::Int, compress::Bool,
                            scalar::Bool, progress::Bool)
    open(filename, "w") do io
        _stream_mzxml_open(io, length(scans))
        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzXML: ") : nothing
        for scan in scans
            _stream_mzxml_spectrum(io, scan;
                                   precision = precision, compress = compress,
                                   scalar = scalar)
            prog === nothing || next!(prog)
        end
        _stream_mzxml_close(io)
    end
    return filename
end


function save_mzxml(filename::AbstractString, scan::MSscans;
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true)
    return _save_mzxml_msscans_vector(filename, [scan];
                                      precision = precision, compress = compress,
                                      scalar = true, progress = progress)
end

function save_mzxml(filename::AbstractString, scans::Vector{MSscans};
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true)
    return _save_mzxml_msscans_vector(filename, scans;
                                      precision = precision, compress = compress,
                                      scalar = false, progress = progress)
end

function _save_mzxml_msscans_vector(filename::AbstractString, scans::Vector{MSscans};
                                    precision::Int, compress::Bool,
                                    scalar::Bool, progress::Bool)
    open(filename, "w") do io
        _stream_mzxml_open(io, length(scans))
        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzXML: ") : nothing
        for sc in scans
            _stream_mzxml_msscans_spectrum(io, sc;
                                           precision = precision, compress = compress,
                                           scalar = scalar)
            prog === nothing || next!(prog)
        end
        _stream_mzxml_close(io)
    end
    return filename
end
