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
* `warn::Bool = true`     — warn when writing a bare scan vector to mzML (see below)

!!! warning "Save the `MSrun`, not a bare vector"
    [`load`](@ref) returns an [`MSrun`](@ref), which carries the file-level
    metadata (instrument, software, source file, data processing) and any
    pre-computed chromatograms. **Slicing or `collect`-ing an `MSrun` yields a
    plain `Vector` and drops that information.** Saving such a bare vector to
    mzML writes only the spectra and emits a warning. Keep (and save) the
    `MSrun` itself to round-trip everything; pass `warn = false` to silence the
    warning when writing bare vectors on purpose.

# Examples
```julia
run = load("input.mzML")                            # an MSrun
save(run, "output.mzML")                            # full round-trip (metadata + chromatograms)
save(run, "output.mzXML"; precision = 32)           # smaller file, lossy
save(run[1], "single_scan.mzML")                    # one scan
save(run, "quiet.mzML"; progress = false)           # no progress bar
save(collect(run), "spectra.mzML")                  # bare vector → warns (metadata dropped)
save(collect(run), "spectra.mzML"; warn = false)    # …silenced
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
# Tee IO that mirrors writes to an underlying file AND a SHA-1 context, so we
# can compute the indexedmzML fileChecksum in a single streaming pass instead
# of re-reading the file at the end.
# ----------------------------------------------------------------------------

mutable struct TeeSHA1IO <: IO
    file::IO
    ctx::SHA.SHA1_CTX
    sealed::Bool   # set true after digest!(ctx); subsequent writes don't hash
    TeeSHA1IO(file::IO) = new(file, SHA.SHA1_CTX(), false)
end

# `write(::IO, args...)` decomposes into single-arg writes, so overriding the
# byte-vector and string cases is enough to cover every call in this module.
# Specific overrides for the concrete byte containers we hit, to avoid the
# `AbstractVector{UInt8}` ↔ `StridedArray` ambiguity in Base.
function Base.write(t::TeeSHA1IO, bytes::Vector{UInt8})
    n = write(t.file, bytes)
    t.sealed || SHA.update!(t.ctx, bytes)
    return n
end
function Base.write(t::TeeSHA1IO, bytes::Base.CodeUnits)
    n = write(t.file, bytes)
    t.sealed || SHA.update!(t.ctx, bytes)
    return n
end
# String / SubString need explicit overrides because Base.write(::IO, ::String)
# takes an unsafe_write fast path that bypasses an AbstractString dispatch.
Base.write(t::TeeSHA1IO, s::String)            = write(t, Vector{UInt8}(codeunits(s)))
Base.write(t::TeeSHA1IO, s::SubString{String}) = write(t, Vector{UInt8}(codeunits(s)))
Base.write(t::TeeSHA1IO, c::Char)              = write(t, string(c))
Base.write(t::TeeSHA1IO, b::UInt8)             = write(t, [b])
Base.position(t::TeeSHA1IO)                    = position(t.file)

# Finalise the hash (call exactly once, right after the `<fileChecksum>` open
# tag has been written). After this call no further bytes are mixed into the
# SHA-1 digest, but writes continue to flow to the underlying file.
function _seal_sha1!(t::TeeSHA1IO)
    @assert !t.sealed "TeeSHA1IO already sealed"
    hash = bytes2hex(SHA.digest!(t.ctx))
    t.sealed = true
    return hash
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


# -- Optional cvParams from MSscan.metadata (mzML round-trip) ----------------
#
# Locations match those parsed by `_read_mzml_extra_metadata` in src/mzml.jl.
# Each emission is gated by `haskey(md, key)` so that absent fields stay
# absent in the output (rather than being filled with zero/empty defaults).

# Direct children of <spectrum>: spectrum_title, lowest/highest observed m/z.
function _stream_mzml_metadata_spectrum_cvs(io::IO, md::Dict{String,Any})
    if haskey(md, "spectrum_title")
        _stream_cvParam(io, CV_SPECTRUM_TITLE, "spectrum title";
                       value = string(md["spectrum_title"]))
    end
    if haskey(md, "lowest_observed_mz")
        _stream_cvParam(io, CV_LOWEST_OBSERVED_MZ, "lowest observed m/z";
                       value = string(md["lowest_observed_mz"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    if haskey(md, "highest_observed_mz")
        _stream_cvParam(io, CV_HIGHEST_OBSERVED_MZ, "highest observed m/z";
                       value = string(md["highest_observed_mz"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
end

# Inside <scanList>/<scan>: mass_resolving_power, ion_injection_time, filter_string.
function _stream_mzml_metadata_scan_cvs(io::IO, md::Dict{String,Any})
    if haskey(md, "mass_resolving_power")
        _stream_cvParam(io, CV_MASS_RESOLVING_POWER, "mass resolving power";
                       value = string(md["mass_resolving_power"]))
    end
    if haskey(md, "ion_injection_time")
        _stream_cvParam(io, CV_ION_INJECTION_TIME, "ion injection time";
                       value = string(md["ion_injection_time"]),
                       unit_cv = "UO", unit_acc = CV_UNIT_MILLISECOND,
                       unit_name = "millisecond")
    end
    if haskey(md, "filter_string")
        _stream_cvParam(io, CV_FILTER_STRING, "filter string";
                       value = string(md["filter_string"]))
    end
end

# Inside <scan>: <scanWindowList> with the lower/upper limits, only when
# either limit is present in metadata.
function _stream_mzml_scan_window(io::IO, md::Dict{String,Any})
    has_lo = haskey(md, "scan_window_lower")
    has_hi = haskey(md, "scan_window_upper")
    (has_lo || has_hi) || return
    write(io, "<scanWindowList count=\"1\">\n<scanWindow>\n")
    if has_lo
        _stream_cvParam(io, CV_SCAN_WINDOW_LOWER, "scan window lower limit";
                       value = string(md["scan_window_lower"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    if has_hi
        _stream_cvParam(io, CV_SCAN_WINDOW_UPPER, "scan window upper limit";
                       value = string(md["scan_window_upper"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    write(io, "</scanWindow>\n</scanWindowList>\n")
end

# Inside <precursor>: <isolationWindow> with target + offsets, when any is
# present in metadata.
function _stream_mzml_isolation_window(io::IO, md::Dict{String,Any})
    has_t = haskey(md, "isolation_window_target_mz")
    has_l = haskey(md, "isolation_window_lower_offset")
    has_u = haskey(md, "isolation_window_upper_offset")
    (has_t || has_l || has_u) || return
    write(io, "<isolationWindow>\n")
    if has_t
        _stream_cvParam(io, CV_ISO_WINDOW_TARGET, "isolation window target m/z";
                       value = string(md["isolation_window_target_mz"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    if has_l
        _stream_cvParam(io, CV_ISO_WINDOW_LOWER, "isolation window lower offset";
                       value = string(md["isolation_window_lower_offset"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    if has_u
        _stream_cvParam(io, CV_ISO_WINDOW_UPPER, "isolation window upper offset";
                       value = string(md["isolation_window_upper_offset"]),
                       unit_cv = "MS", unit_acc = CV_UNIT_MZ, unit_name = "m/z")
    end
    write(io, "</isolationWindow>\n")
end

# Inside <selectedIon>: peak intensity of the selected precursor.
function _stream_mzml_selected_ion_intensity(io::IO, md::Dict{String,Any})
    haskey(md, "selected_ion_peak_intensity") || return
    _stream_cvParam(io, CV_SELECTED_PEAK_INT, "peak intensity";
                   value = string(md["selected_ion_peak_intensity"]),
                   unit_cv = "MS", unit_acc = "MS:1000131",
                   unit_name = "number of detector counts")
end


# ----------------------------------------------------------------------------
# Streaming mzML writer
# ----------------------------------------------------------------------------

function _stream_mzml_open(io::IO, scancount::Int;
                           file_metadata::Dict{String,Any} = Dict{String,Any}(),
                           with_xml_decl::Bool = true)
    if with_xml_decl
        write(io, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
    end
    write(io, "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"",
          " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
          " version=\"1.1.0\">\n")

    # cvList is fixed (MS + UO) — same in any mzML file we emit.
    write(io, "<cvList count=\"2\">\n",
              "<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\"",
              " URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"",
              " version=\"4.1.0\"/>\n",
              "<cv id=\"UO\" fullName=\"Unit Ontology\"",
              " URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"",
              " version=\"12:10:2011\"/>\n",
              "</cvList>\n")

    # File description, source files, param groups, software, instruments,
    # data processing — emitted from parsed metadata when available, falls
    # back to a minimal stub otherwise.
    _stream_mzml_file_metadata(io, file_metadata)

    # Pick the first instrument id (if any) and first data-processing id for
    # the run's defaultInstrumentConfigurationRef / defaultDataProcessingRef.
    instr_id = _first_id(get(file_metadata, "instruments", nothing), "IC1")
    dp_id    = _first_id(get(file_metadata, "data_processing", nothing), "MassJExport")
    write(io, "<run id=\"run1\" defaultInstrumentConfigurationRef=\"", instr_id, "\">\n",
              "<spectrumList count=\"", string(scancount),
              "\" defaultDataProcessingRef=\"", dp_id, "\">\n")
end

# Pick the id of the first dict in a vector, with a fallback when the vector
# is missing or empty. Used to wire up the run element's default refs.
function _first_id(v, fallback::String)
    v === nothing && return fallback
    isempty(v) && return fallback
    id = get(v[1], "id", nothing)
    id === nothing ? fallback : id
end


function _stream_mzml_close(io::IO;
                            chromatograms::Vector{Chromatogram} = Chromatogram[],
                            precision::Int = 64, compress::Bool = true)
    write(io, "</spectrumList>\n")
    _stream_mzml_chromatogramList(io, chromatograms;
                                  precision = precision, compress = compress)
    write(io, "</run>\n</mzML>\n")
end


function _stream_mzml_spectrum(io::IO, scan::MSscan, index::Int;
                               precision::Int = 64, compress::Bool = true,
                               scalar::Bool = false)
    write(io, "<spectrum index=\"", string(index),
              "\" id=\"scan=", string(scan.num),
              "\" defaultArrayLength=\"", string(length(scan.mz)), "\">\n")

    # mzML 1.1 schema requires: cvParam* userParam* (then child elements).
    # All cvParams come first; any userParams emitted strictly after.
    _stream_cvParam(io, CV_MS_LEVEL, "ms level"; value = string(scan.level))
    if scan.level >= 2
        _stream_cvParam(io, CV_MSN_SPECTRUM, "MSn spectrum")
    end

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
    _stream_mzml_metadata_spectrum_cvs(io, scan.metadata)

    # userParams after cvParams (schema order).
    if scalar
        _stream_userParam(io, MASSJ_SCALAR_PARAM; value = "true")
    end

    write(io, "<scanList count=\"1\">\n")
    _stream_cvParam(io, "MS:1000795", "no combination")
    write(io, "<scan>\n")
    _stream_cvParam(io, CV_SCAN_START_TIME, "scan start time";
                   value = string(scan.rt), unit_cv = "UO",
                   unit_acc = CV_UNIT_MINUTE, unit_name = "minute")
    _stream_mzml_metadata_scan_cvs(io, scan.metadata)
    _stream_mzml_scan_window(io, scan.metadata)
    write(io, "</scan>\n</scanList>\n")

    if scan.level >= 2 && scan.precursor > 0
        write(io, "<precursorList count=\"1\">\n<precursor>\n")
        _stream_mzml_isolation_window(io, scan.metadata)
        write(io, "<selectedIonList count=\"1\">\n<selectedIon>\n")
        _stream_cvParam(io, CV_SELECTED_ION_MZ, "selected ion m/z";
                       value = string(scan.precursor),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        if scan.chargeState != 0
            _stream_cvParam(io, CV_CHARGE_STATE, "charge state";
                           value = string(scan.chargeState))
        end
        _stream_mzml_selected_ion_intensity(io, scan.metadata)
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
    # The original first scan number is preserved in MassJ:num (see below);
    # the spectrum's `id` attribute is derived from the index so the schema's
    # uniqueness constraint on spectrum ids holds when several MSscans are
    # saved together.
    write(io, "<spectrum index=\"", string(index),
              "\" id=\"scan=", string(index + 1),
              "\" defaultArrayLength=\"", string(length(scan.mz)), "\">\n")

    # mzML 1.1 schema requires: cvParam* userParam* (then child elements).
    # All cvParams come first; userParams are emitted after the metadata cvs.
    lvl = isempty(scan.level) ? 1 : scan.level[1]
    _stream_cvParam(io, CV_MS_LEVEL, "ms level"; value = string(lvl))
    if lvl >= 2
        _stream_cvParam(io, CV_MSN_SPECTRUM, "MSn spectrum")
    end

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
    _stream_mzml_metadata_spectrum_cvs(io, scan.metadata)

    # userParams (schema-mandated to come after cvParams).
    _stream_userParam(io, MASSJ_CONTAINER_PARAM; value = "MSscans")
    if scalar
        _stream_userParam(io, MASSJ_SCALAR_PARAM; value = "true")
    end
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

    rt0 = isempty(scan.rt) ? 0.0 : scan.rt[1]
    write(io, "<scanList count=\"1\">\n")
    _stream_cvParam(io, "MS:1000795", "no combination")
    write(io, "<scan>\n")
    _stream_cvParam(io, CV_SCAN_START_TIME, "scan start time";
                   value = string(rt0), unit_cv = "UO",
                   unit_acc = CV_UNIT_MINUTE, unit_name = "minute")
    _stream_mzml_metadata_scan_cvs(io, scan.metadata)
    _stream_mzml_scan_window(io, scan.metadata)
    write(io, "</scan>\n</scanList>\n")

    prec0 = isempty(scan.precursor) ? 0.0 : scan.precursor[1]
    if lvl >= 2 && prec0 > 0
        write(io, "<precursorList count=\"1\">\n<precursor>\n")
        _stream_mzml_isolation_window(io, scan.metadata)
        write(io, "<selectedIonList count=\"1\">\n<selectedIon>\n")
        _stream_cvParam(io, CV_SELECTED_ION_MZ, "selected ion m/z";
                       value = string(prec0),
                       unit_cv = "MS", unit_acc = "MS:1000040", unit_name = "m/z")
        chg0 = isempty(scan.chargeState) ? 0 : scan.chargeState[1]
        if chg0 != 0
            _stream_cvParam(io, CV_CHARGE_STATE, "charge state"; value = string(chg0))
        end
        _stream_mzml_selected_ion_intensity(io, scan.metadata)
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

# Warn when a bare scan vector (rather than the `MSrun` returned by `load`) is
# written to mzML. Run-level metadata and pre-computed chromatograms live on the
# `MSrun` wrapper, not on a scan vector, so they are silently absent from the
# output; slicing or `collect`-ing an `MSrun` produces such a bare vector. Fires
# on every lossy save (deterministic, easy to test); silence with `warn = false`.
function _warn_dropped_run_metadata(n::Integer; warn::Bool)
    warn || return nothing
    @warn "save: writing a bare vector of $n spectra to mzML. Run-level " *
          "metadata (instrument, software, source file, data processing) and " *
          "any pre-computed chromatograms are NOT written: they live on the " *
          "`MSrun` returned by `load`, not on a scan vector, and are dropped " *
          "when an `MSrun` is sliced or `collect`-ed. To preserve them, save " *
          "the original `MSrun`, or — after processing its scans — wrap them " *
          "in an `MSrun` that carries the original's metadata: " *
          "`save(MSrun(scans, run), file)` (equivalently " *
          "`MSrun(scans, run.metadata, run.chromatograms)`). Pass " *
          "`warn = false` to silence this message."
    return nothing
end

"""
    save_mzml(filename::AbstractString, data;
              precision::Int = 64, compress::Bool = true,
              progress::Bool = true, warn::Bool = true) -> filename
Write a [`MSscan`](@ref), [`MSscans`](@ref), `Vector{MSscan}` /
`Vector{MSscans}`, or an [`MSrun`](@ref) to an mzML file. The file is written
one spectrum at a time; peak RAM is bounded by the largest single spectrum, not
the total file size.

When `progress = true` (default) a `ProgressMeter` bar is shown while writing.

!!! warning "Run-level metadata"
    File-level metadata (instrument, software, source file, data processing) and
    pre-computed chromatograms are carried by the [`MSrun`](@ref) that [`load`](@ref)
    returns, **not** by a plain scan vector. Saving a bare `Vector{MSscan}` /
    `Vector{MSscans}` therefore drops that information and emits a warning; save
    the `MSrun` itself to preserve it. Set `warn = false` to silence the warning
    (e.g. in batch scripts that intentionally write bare vectors).

Round-trips through [`load`](@ref): the loaded value has the same type as the
saved one (scalar `MSscan` / `MSscans` come back bare; vectors come back as
vectors; an `MSrun` comes back as an `MSrun`).
"""
function save_mzml(filename::AbstractString, scans::Vector{MSscan};
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true, indexed::Bool = true,
                   warn::Bool = true)
    _warn_dropped_run_metadata(length(scans); warn = warn)
    return _save_mzml_vector(filename, scans;
                             precision = precision, compress = compress,
                             scalar = false, progress = progress,
                             file_metadata = Dict{String,Any}(),
                             indexed = indexed)
end

# MSrun preserves file-level metadata + chromatograms across the round-trip.
function save_mzml(filename::AbstractString, run::MSrun;
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true, indexed::Bool = true,
                   warn::Bool = true)
    return _save_mzml_vector(filename, run.scans;
                             precision = precision, compress = compress,
                             scalar = false, progress = progress,
                             file_metadata = run.metadata,
                             chromatograms = run.chromatograms,
                             indexed = indexed)
end

function save_mzml(filename::AbstractString, scan::MSscan;
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true, indexed::Bool = true,
                   warn::Bool = true)
    return _save_mzml_vector(filename, [scan];
                             precision = precision, compress = compress,
                             scalar = true, progress = progress,
                             file_metadata = Dict{String,Any}(),
                             indexed = indexed)
end

function _save_mzml_vector(filename::AbstractString, scans::Vector{MSscan};
                           precision::Int, compress::Bool,
                           scalar::Bool, progress::Bool,
                           file_metadata::Dict{String,Any} = Dict{String,Any}(),
                           chromatograms::Vector{Chromatogram} = Chromatogram[],
                           indexed::Bool = true)
    open(filename, "w") do raw_io
        io = indexed ? TeeSHA1IO(raw_io) : raw_io

        if indexed
            write(io, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
                      "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\"",
                      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
                      " xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml",
                      " http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd\">\n")
            _stream_mzml_open(io, length(scans);
                              file_metadata = file_metadata, with_xml_decl = false)
        else
            _stream_mzml_open(io, length(scans); file_metadata = file_metadata)
        end

        spec_offsets  = Vector{Tuple{String,Int}}(undef, 0)
        chrom_offsets = Vector{Tuple{String,Int}}(undef, 0)

        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzML: ") : nothing
        for (i, scan) in enumerate(scans)
            indexed && push!(spec_offsets, ("scan=$(scan.num)", position(raw_io)))
            _stream_mzml_spectrum(io, scan, i - 1;
                                  precision = precision, compress = compress,
                                  scalar = scalar)
            prog === nothing || next!(prog)
        end
        write(io, "</spectrumList>\n")

        # chromatogramList — inlined with per-chromatogram offset capture so
        # the indexed wrapper can list them too.
        if !isempty(chromatograms)
            write(io, "<chromatogramList count=\"", string(length(chromatograms)),
                      "\" defaultDataProcessingRef=\"MassJExport\">\n")
            for (i, c) in enumerate(chromatograms)
                id = "chrom_" * string(i)
                indexed && push!(chrom_offsets, (id, position(raw_io)))
                write(io, "<chromatogram id=\"", id,
                          "\" index=\"", string(i - 1),
                          "\" defaultArrayLength=\"", string(length(c.rt)), "\">\n")
                _stream_cvParam(io, CV_TIC_CHROMATOGRAM, "total ion current chromatogram")
                write(io, "<binaryDataArrayList count=\"2\">\n")
                _stream_mzml_chrom_binary(io, c.rt,  :time;
                                          precision = precision, compress = compress)
                _stream_mzml_chrom_binary(io, c.ic,  :intensity;
                                          precision = precision, compress = compress)
                write(io, "</binaryDataArrayList>\n</chromatogram>\n")
            end
            write(io, "</chromatogramList>\n")
        end

        write(io, "</run>\n</mzML>\n")

        if indexed
            _stream_indexed_footer(io, raw_io, spec_offsets, chrom_offsets)
        end
    end
    return filename
end


function save_mzml(filename::AbstractString, scan::MSscans;
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true, indexed::Bool = true,
                   warn::Bool = true)
    return _save_mzml_msscans_vector(filename, [scan];
                                     precision = precision, compress = compress,
                                     scalar = true, progress = progress,
                                     indexed = indexed)
end

function save_mzml(filename::AbstractString, scans::Vector{MSscans};
                   precision::Int = 64, compress::Bool = true,
                   progress::Bool = true, indexed::Bool = true,
                   warn::Bool = true)
    _warn_dropped_run_metadata(length(scans); warn = warn)
    return _save_mzml_msscans_vector(filename, scans;
                                     precision = precision, compress = compress,
                                     scalar = false, progress = progress,
                                     indexed = indexed)
end

function _save_mzml_msscans_vector(filename::AbstractString, scans::Vector{MSscans};
                                   precision::Int, compress::Bool,
                                   scalar::Bool, progress::Bool,
                                   file_metadata::Dict{String,Any} = Dict{String,Any}(),
                                   indexed::Bool = true)
    open(filename, "w") do raw_io
        io = indexed ? TeeSHA1IO(raw_io) : raw_io

        if indexed
            write(io, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
                      "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\"",
                      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
                      " xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml",
                      " http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd\">\n")
            _stream_mzml_open(io, length(scans);
                              file_metadata = file_metadata, with_xml_decl = false)
        else
            _stream_mzml_open(io, length(scans); file_metadata = file_metadata)
        end

        spec_offsets = Vector{Tuple{String,Int}}(undef, 0)

        prog = progress && length(scans) > 1 ?
            Progress(length(scans); desc = "Writing mzML: ") : nothing
        for (i, sc) in enumerate(scans)
            # Spectrum id matches what _stream_mzml_msscans_spectrum writes
            # (scan=<index+1>, see the uniqueness comment there).
            indexed && push!(spec_offsets, ("scan=" * string(i), position(raw_io)))
            _stream_mzml_msscans_spectrum(io, sc, i - 1;
                                          precision = precision, compress = compress,
                                          scalar = scalar)
            prog === nothing || next!(prog)
        end
        write(io, "</spectrumList>\n</run>\n</mzML>\n")

        if indexed
            _stream_indexed_footer(io, raw_io, spec_offsets,
                                   Vector{Tuple{String,Int}}())
        end
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
Write a [`MSscan`](@ref), [`MSscans`](@ref), `Vector{MSscan}` /
`Vector{MSscans}`, or an [`MSrun`](@ref) to an mzXML file. Streams one spectrum
at a time; peak RAM is bounded by the largest single spectrum. mzXML interleaves
m/z and intensity in a single `<peaks>` blob and uses big-endian (network) byte
order.

!!! note
    The mzXML writer emits spectra only; file-level metadata and chromatograms
    are not part of the output regardless of the input type. An [`MSrun`](@ref)
    is accepted (its scans are written) so that code keeping the `MSrun` from
    [`load`](@ref) — the recommended pattern for mzML — also works here. The
    `warn` keyword is accepted for signature parity with [`save_mzml`](@ref) but
    has no effect.
"""
function save_mzxml(filename::AbstractString, scans::Vector{MSscan};
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true, warn::Bool = true)
    return _save_mzxml_vector(filename, scans;
                              precision = precision, compress = compress,
                              scalar = false, progress = progress)
end

# Accepted so the MSrun returned by `load` can be saved as mzXML without a
# MethodError; only its scans are written (mzXML carries no run-level metadata).
function save_mzxml(filename::AbstractString, run::MSrun;
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true, warn::Bool = true)
    return _save_mzxml_vector(filename, run.scans;
                              precision = precision, compress = compress,
                              scalar = false, progress = progress)
end

function save_mzxml(filename::AbstractString, scan::MSscan;
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true, warn::Bool = true)
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
                    progress::Bool = true, warn::Bool = true)
    return _save_mzxml_msscans_vector(filename, [scan];
                                      precision = precision, compress = compress,
                                      scalar = true, progress = progress)
end

function save_mzxml(filename::AbstractString, scans::Vector{MSscans};
                    precision::Int = 64, compress::Bool = true,
                    progress::Bool = true, warn::Bool = true)
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


# ============================================================================
# File-level metadata writing (matches _read_mzml_file_metadata)
# ============================================================================

# Emit a cvParam from a parsed Dict (the shape produced by _cvparams_to_dicts).
function _stream_cvParam_from_dict(io::IO, d::AbstractDict)
    acc  = get(d, "accession", "")
    nm   = get(d, "name", "")
    isempty(acc) && return
    _stream_cvParam(io, acc, nm;
                   value     = get(d, "value", ""),
                   unit_cv   = get(d, "unitCvRef", ""),
                   unit_acc  = get(d, "unitAccession", ""),
                   unit_name = get(d, "unitName", ""))
end

# Top-level: <fileDescription>, <referenceableParamGroupList>, <softwareList>,
# <instrumentConfigurationList>, <dataProcessingList>. Falls back to the
# original minimal stub when `file_metadata` is empty.
function _stream_mzml_file_metadata(io::IO, md::Dict{String,Any})
    isempty(md) && return _stream_mzml_minimal_stub(io)

    # -- fileDescription ----------------------------------------------------
    write(io, "<fileDescription>\n<fileContent>\n")
    fc = get(md, "file_content", nothing)
    if fc !== nothing && !isempty(fc)
        for cv in fc
            _stream_cvParam_from_dict(io, cv)
        end
    else
        _stream_cvParam(io, "MS:1000579", "MS1 spectrum")
    end
    write(io, "</fileContent>\n")

    sfs = get(md, "source_files", nothing)
    if sfs !== nothing && !isempty(sfs)
        write(io, "<sourceFileList count=\"", string(length(sfs)), "\">\n")
        for sf in sfs
            write(io, "<sourceFile")
            for attr in ("id", "name", "location")
                v = get(sf, attr, "")
                isempty(v) || write(io, " ", attr, "=\"", _xmlescape(v), "\"")
            end
            write(io, ">\n")
            for cv in get(sf, "cv_params", [])
                _stream_cvParam_from_dict(io, cv)
            end
            write(io, "</sourceFile>\n")
        end
        write(io, "</sourceFileList>\n")
    end
    write(io, "</fileDescription>\n")

    # -- referenceableParamGroupList ----------------------------------------
    pgs = get(md, "param_groups", nothing)
    if pgs !== nothing && !isempty(pgs)
        write(io, "<referenceableParamGroupList count=\"", string(length(pgs)), "\">\n")
        for pg in pgs
            id = get(pg, "id", "")
            write(io, "<referenceableParamGroup id=\"", _xmlescape(id), "\">\n")
            for cv in get(pg, "cv_params", [])
                _stream_cvParam_from_dict(io, cv)
            end
            write(io, "</referenceableParamGroup>\n")
        end
        write(io, "</referenceableParamGroupList>\n")
    end

    # -- softwareList -------------------------------------------------------
    # Always include a MassJ software entry so any data_processing stub
    # referencing softwareRef="MassJ" satisfies the schema's keyref.
    sws = get(md, "software", nothing)
    has_massj = sws !== nothing &&
                any(sw -> get(sw, "id", "") == "MassJ", sws)
    base_sws = sws === nothing ? Dict{String,Any}[] : sws
    n_total  = length(base_sws) + (has_massj ? 0 : 1)

    write(io, "<softwareList count=\"", string(n_total), "\">\n")
    for sw in base_sws
        write(io, "<software")
        for attr in ("id", "version")
            v = get(sw, attr, "")
            isempty(v) || write(io, " ", attr, "=\"", _xmlescape(v), "\"")
        end
        write(io, ">\n")
        for cv in get(sw, "cv_params", [])
            _stream_cvParam_from_dict(io, cv)
        end
        write(io, "</software>\n")
    end
    if !has_massj
        write(io, "<software id=\"MassJ\" version=\"0.1\">\n")
        _stream_cvParam(io, "MS:1000799", "custom unreleased software tool"; value = "MassJ")
        write(io, "</software>\n")
    end
    write(io, "</softwareList>\n")

    # -- instrumentConfigurationList ----------------------------------------
    instrs = get(md, "instruments", nothing)
    if instrs !== nothing && !isempty(instrs)
        write(io, "<instrumentConfigurationList count=\"", string(length(instrs)), "\">\n")
        for ic in instrs
            id = get(ic, "id", "")
            write(io, "<instrumentConfiguration id=\"", _xmlescape(id), "\">\n")
            for ref in get(ic, "param_group_refs", String[])
                write(io, "<referenceableParamGroupRef ref=\"", _xmlescape(ref), "\"/>\n")
            end
            for cv in get(ic, "cv_params", [])
                _stream_cvParam_from_dict(io, cv)
            end
            comps = get(ic, "components", [])
            if !isempty(comps)
                write(io, "<componentList count=\"", string(length(comps)), "\">\n")
                for c in comps
                    t   = get(c, "type", "source")
                    ord = get(c, "order", "1")
                    write(io, "<", t, " order=\"", _xmlescape(ord), "\">\n")
                    for cv in get(c, "cv_params", [])
                        _stream_cvParam_from_dict(io, cv)
                    end
                    write(io, "</", t, ">\n")
                end
                write(io, "</componentList>\n")
            end
            swRef = get(ic, "software_ref", "")
            isempty(swRef) || write(io, "<softwareRef ref=\"", _xmlescape(swRef), "\"/>\n")
            write(io, "</instrumentConfiguration>\n")
        end
        write(io, "</instrumentConfigurationList>\n")
    else
        write(io, "<instrumentConfigurationList count=\"1\">\n",
                  "<instrumentConfiguration id=\"IC1\">\n")
        _stream_cvParam(io, "MS:1000031", "instrument model")
        write(io, "</instrumentConfiguration>\n</instrumentConfigurationList>\n")
    end

    # -- dataProcessingList -------------------------------------------------
    dps = get(md, "data_processing", nothing)
    if dps !== nothing && !isempty(dps)
        write(io, "<dataProcessingList count=\"", string(length(dps)), "\">\n")
        for dp in dps
            id = get(dp, "id", "")
            write(io, "<dataProcessing id=\"", _xmlescape(id), "\">\n")
            for m in get(dp, "methods", [])
                ord = get(m, "order", "0")
                swR = get(m, "softwareRef", "")
                write(io, "<processingMethod order=\"", _xmlescape(ord), "\"")
                isempty(swR) || write(io, " softwareRef=\"", _xmlescape(swR), "\"")
                write(io, ">\n")
                for cv in get(m, "cv_params", [])
                    _stream_cvParam_from_dict(io, cv)
                end
                write(io, "</processingMethod>\n")
            end
            write(io, "</dataProcessing>\n")
        end
        write(io, "</dataProcessingList>\n")
    else
        write(io, "<dataProcessingList count=\"1\">\n",
                  "<dataProcessing id=\"MassJExport\">\n",
                  "<processingMethod order=\"0\" softwareRef=\"MassJ\">\n")
        _stream_cvParam(io, "MS:1000544", "Conversion to mzML")
        write(io, "</processingMethod>\n</dataProcessing>\n</dataProcessingList>\n")
    end
end

# Minimal placeholder block used when no metadata is supplied (Vector{MSscan}
# save without an MSrun). Identical to the pre-MSrun output.
function _stream_mzml_minimal_stub(io::IO)
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
end


# ----------------------------------------------------------------------------
# Chromatogram block emission (Phase 3)
# ----------------------------------------------------------------------------

# Emit a <chromatogramList> after </spectrumList> when run.chromatograms is
# non-empty. Each Chromatogram is written as a TIC chromatogram with two
# binary arrays: time (UO:0000031, minutes) and intensity.
function _stream_mzml_chromatogramList(io::IO, chroms::Vector{Chromatogram};
                                       precision::Int = 64, compress::Bool = true)
    isempty(chroms) && return
    write(io, "<chromatogramList count=\"", string(length(chroms)), "\""
              * " defaultDataProcessingRef=\"MassJExport\">\n")
    for (i, c) in enumerate(chroms)
        write(io, "<chromatogram id=\"chrom_", string(i),
                  "\" index=\"", string(i - 1),
                  "\" defaultArrayLength=\"", string(length(c.rt)), "\">\n")
        _stream_cvParam(io, CV_TIC_CHROMATOGRAM, "total ion current chromatogram")
        write(io, "<binaryDataArrayList count=\"2\">\n")
        _stream_mzml_chrom_binary(io, c.rt,  :time;
                                  precision = precision, compress = compress)
        _stream_mzml_chrom_binary(io, c.ic,  :intensity;
                                  precision = precision, compress = compress)
        write(io, "</binaryDataArrayList>\n</chromatogram>\n")
    end
    write(io, "</chromatogramList>\n")
end

function _stream_mzml_chrom_binary(io::IO, data::Vector{Float64}, kind::Symbol;
                                   precision::Int = 64, compress::Bool = true)
    b64, _ = _encode_binary(data; precision = precision, compress = compress,
                            endian = :little)
    write(io, "<binaryDataArray encodedLength=\"", string(length(b64)), "\">\n")
    _stream_cvParam(io, precision == 64 ? CV_64BIT : CV_32BIT,
                   precision == 64 ? "64-bit float" : "32-bit float")
    _stream_cvParam(io, compress ? CV_ZLIB : CV_NO_COMPRESSION,
                   compress ? "zlib compression" : "no compression")
    if kind === :time
        _stream_cvParam(io, CV_TIME_ARRAY, "time array";
                       unit_cv = "UO", unit_acc = CV_UNIT_MINUTE,
                       unit_name = "minute")
    else
        _stream_cvParam(io, CV_INT_ARRAY, "intensity array";
                       unit_cv = "MS", unit_acc = "MS:1000131",
                       unit_name = "number of detector counts")
    end
    write(io, "<binary>", b64, "</binary>\n</binaryDataArray>\n")
end


# ----------------------------------------------------------------------------
# Indexed-mzML footer: <indexList>, <indexListOffset>, <fileChecksum>
# ----------------------------------------------------------------------------

"""
    _stream_indexed_footer(tee, raw, spec_offsets, chrom_offsets)
Append the `<indexList>`, `<indexListOffset>`, and `<fileChecksum>` blocks
that complete an indexed-mzML file. Per the mzML 1.1 indexed schema:

* the SHA-1 in `<fileChecksum>` is taken over every byte from the start of
  the file up to **and including** the `<fileChecksum>` open tag — which is
  why we tee the file through a [`TeeSHA1IO`](@ref) and seal the digest
  immediately after writing that open tag.
* `<indexListOffset>` is the byte offset (in the raw file) where the
  `<indexList>` element starts.
"""
function _stream_indexed_footer(io::TeeSHA1IO, raw_io::IO,
                                spec_offsets::Vector{Tuple{String,Int}},
                                chrom_offsets::Vector{Tuple{String,Int}})
    # Capture <indexList> start BEFORE writing it.
    indexlist_offset = position(raw_io)

    n_indices = 1 + (isempty(chrom_offsets) ? 0 : 1)
    write(io, "<indexList count=\"", string(n_indices), "\">\n")

    write(io, "<index name=\"spectrum\">\n")
    for (id, off) in spec_offsets
        write(io, "<offset idRef=\"", _xmlescape(id), "\">", string(off), "</offset>\n")
    end
    write(io, "</index>\n")

    if !isempty(chrom_offsets)
        write(io, "<index name=\"chromatogram\">\n")
        for (id, off) in chrom_offsets
            write(io, "<offset idRef=\"", _xmlescape(id), "\">", string(off), "</offset>\n")
        end
        write(io, "</index>\n")
    end

    write(io, "</indexList>\n",
              "<indexListOffset>", string(indexlist_offset), "</indexListOffset>\n",
              "<fileChecksum>")

    # SHA-1 is sealed *here* — after the <fileChecksum> open tag is written
    # (still part of the hashed content) but before its text content.
    hash_hex = _seal_sha1!(io)

    # Remaining bytes are written to the underlying file only; the digest
    # is finalised so SHA.update! is no longer called.
    write(io, hash_hex, "</fileChecksum>\n</indexedmzML>\n")
end
