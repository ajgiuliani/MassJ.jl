# =============================================================================
# Fast DOM-free streaming reader for imzML.
#
# The approach used here — parsing the bulk <spectrumList> as a text stream
# instead of building a full XML DOM, reading binary arrays directly from the
# companion .ibd file — is adapted from julia_mzML_imzML by Dr Robert Winkler:
#
#     https://codeberg.org/LabABI/julia_mzML_imzML
#     Copyright (c) 2023 Dr Robert Winkler — MIT License
#     (full text in THIRD_PARTY_NOTICES.md at the repository root)
#
# MassJ reimplements the technique with a line-oriented parser (robust to
# variable-width values) that populates full MSscans structs, and cites:
#
#     I. Rosas-Román, H. Guillén-Alonso, A. Moreno-Pedraza, R. Winkler.
#     Anal. Chem. 2024. DOI: 10.1021/acs.analchem.3c05853
#
# MassJ as a whole is distributed under GPL-3.0-or-later; this file incorporates
# MIT-licensed ideas whose notice is retained above and in THIRD_PARTY_NOTICES.md.
# =============================================================================


# Small attribute extractors. imzML written by ProteoWizard/etc. puts one
# cvParam per line, so a per-line regex is both sufficient and fast.
const _IMZML_RE_VALUE = r"value=\"([^\"]*)\""
const _IMZML_RE_ACC   = r"accession=\"([^\"]+)\""
const _IMZML_RE_REF   = r"ref=\"([^\"]+)\""
const _IMZML_RE_UNIT  = r"unitAccession=\"([^\"]+)\""
const _IMZML_RE_COUNT = r"count=\"(\d+)\""

@inline function _imzml_attr(line::AbstractString, re::Regex)
    m = match(re, line)
    m === nothing ? nothing : m.captures[1]
end


"""
    load_imzml_stream(filename::String; progress::Bool=true) -> MSrun
Fast DOM-free imzML reader. Parses only the (small) header with LightXML to
resolve per-file invariants — storage mode and the array-type / precision /
compression flags of each `referenceableParamGroup` — then streams the
`<spectrumList>` line by line, reading binary arrays from the `.ibd` file.

Produces the same full `MSscans` structs as [`load_imzml_all_dom`](@ref).
Throws [`ImzmlStreamFallback`](@ref) for layouts it does not handle (no
`spectrumList`, inline binary, precursor data, spectrum count mismatch, …) so
that [`load_imzml_all`](@ref) can fall back to the DOM reader.
"""
function load_imzml_stream(filename::String; progress::Bool=true)
    # Locate the companion .ibd (same resolution as the DOM reader)
    basepath = filename[1:findlast('.', filename)-1]
    ibd_path = basepath * ".ibd"
    if !isfile(ibd_path)
        for ext in (".ibd", ".IBD", ".Ibd")
            candidate = basepath * ext
            if isfile(candidate)
                ibd_path = candidate
                break
            end
        end
    end
    isfile(ibd_path) || throw(ImzmlStreamFallback("companion .ibd file not found for $filename"))

    io = open(filename, "r")
    try
        # ---- Header: accumulate lines up to and including <spectrumList ...> ----
        header = IOBuffer()
        scanCount = 0
        found = false
        while !eof(io)
            line = readline(io; keep=true)
            write(header, line)
            if occursin("<spectrumList", line)
                c = _imzml_attr(line, _IMZML_RE_COUNT)
                c === nothing && throw(ImzmlStreamFallback("<spectrumList> has no count attribute"))
                scanCount = parse(Int, c)
                found = true
                break
            end
        end
        found || throw(ImzmlStreamFallback("no <spectrumList> element found"))

        # Parse just the header by closing the truncated fragment into valid XML.
        header_str = String(take!(header)) * "\n</spectrumList></run></mzML>"
        hdoc = parse_string(header_str)
        local is_continuous::Bool
        local group_flags
        try
            hroot = find_mzml_root(hdoc)
            group_flags = precompute_group_flags(parse_referenceable_param_groups(hroot))
            is_continuous = false
            fileDesc = find_element(hroot, "fileDescription")
            if fileDesc !== nothing
                fileContent = find_element(fileDesc, "fileContent")
                if fileContent !== nothing && has_cv_param(fileContent, CV_IMS_CONTINUOUS)
                    is_continuous = true
                end
            end
        finally
            free(hdoc)
        end

        # ---- Stream the spectrumList ----
        scans = Vector{MSscans}(undef, scanCount)
        prog = (progress && scanCount > 1) ? Progress(scanCount; desc = "Loading imzML: ") : nothing

        open(ibd_path, "r") do ibd_io
            skip(ibd_io, 16)                  # 16-byte UUID at the start of the .ibd
            bytebuf = UInt8[]
            shared_mz = Float64[]
            shared_mz_read = false
            idx = 0

            while !eof(io)
                line = readline(io)
                (occursin("<spectrum ", line) || occursin("<spectrum>", line)) || continue

                idx += 1
                idx <= scanCount || throw(ImzmlStreamFallback("more spectra than declared count ($scanCount)"))

                # Per-spectrum fields (defaults match the DOM reader)
                msLevel = 1
                spectrumType = :unknown
                polarity = ""
                tic = 0.0; basePeakMz = 0.0; basePeakIntensity = 0.0
                rt = 0.0; pos_x = 0; pos_y = 0; pos_z = 0
                mz = Float64[]; int_arr = Float64[]

                # Current <binaryDataArray> accumulator
                cur_active = false
                cur_is_mz = false; cur_is_int = false; cur_ext = false
                cur_64 = false; cur_zlib = false
                cur_off = 0; cur_len = 0; cur_enc = 0

                while true
                    eof(io) && throw(ImzmlStreamFallback("unexpected EOF inside a <spectrum>"))
                    l = readline(io)

                    if occursin("</spectrum>", l)
                        break
                    elseif occursin("<precursorList", l)
                        # MS/MS imaging: defer to the general DOM reader
                        throw(ImzmlStreamFallback("precursor data present"))
                    elseif occursin("<binaryDataArray", l) && !occursin("<binaryDataArrayList", l)
                        cur_active = true
                        cur_is_mz = false; cur_is_int = false; cur_ext = false
                        cur_64 = false; cur_zlib = false
                        cur_off = 0; cur_len = 0; cur_enc = 0
                    elseif occursin("</binaryDataArray>", l)
                        if (cur_is_mz || cur_is_int) && cur_ext && cur_len > 0
                            if is_continuous && cur_is_mz && shared_mz_read
                                mz = shared_mz                 # share the m/z axis by reference
                            else
                                out = read_ibd_array!(ibd_io, cur_off, cur_enc, cur_len,
                                                      cur_64, cur_zlib, bytebuf)
                                cur_is_mz ? (mz = out) : (int_arr = out)
                            end
                        end
                        cur_active = false
                    elseif occursin("referenceableParamGroupRef", l)
                        ref = _imzml_attr(l, _IMZML_RE_REF)
                        if ref !== nothing
                            g = get(group_flags, ref, nothing)
                            if g !== nothing
                                cur_is_mz |= g.is_mz; cur_is_int |= g.is_int; cur_ext |= g.is_external
                                cur_64 |= g.is_64bit; cur_zlib |= g.is_zlib
                            end
                        end
                    elseif occursin("cvParam", l)
                        acc = _imzml_attr(l, _IMZML_RE_ACC)
                        acc === nothing && continue
                        if cur_active
                            if acc == CV_MZ_ARRAY
                                cur_is_mz = true
                            elseif acc == CV_INT_ARRAY
                                cur_is_int = true
                            elseif acc == CV_64BIT
                                cur_64 = true
                            elseif acc == CV_ZLIB
                                cur_zlib = true
                            elseif acc == CV_IMS_EXTERNAL_DATA
                                cur_ext = true
                            elseif acc == CV_IMS_EXT_OFFSET
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (cur_off = parse(Int, v))
                            elseif acc == CV_IMS_EXT_LENGTH
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (cur_len = parse(Int, v))
                            elseif acc == CV_IMS_EXT_ENC_LEN
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (cur_enc = parse(Int, v))
                            end
                        else
                            if acc == CV_MS_LEVEL
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (msLevel = parse(Int, v))
                            elseif acc == CV_CENTROID
                                spectrumType = :centroid
                            elseif acc == CV_PROFILE
                                spectrumType = :profile
                            elseif acc == CV_POSITIVE_SCAN
                                polarity = "+"
                            elseif acc == CV_NEGATIVE_SCAN
                                polarity = "-"
                            elseif acc == CV_TIC
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (tic = parse(Float64, v))
                            elseif acc == CV_BASE_PEAK_MZ
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (basePeakMz = parse(Float64, v))
                            elseif acc == CV_BASE_PEAK_INT
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (basePeakIntensity = parse(Float64, v))
                            elseif acc == CV_SCAN_START_TIME
                                v = _imzml_attr(l, _IMZML_RE_VALUE)
                                if v !== nothing
                                    rtv = parse(Float64, v)
                                    _imzml_attr(l, _IMZML_RE_UNIT) == CV_UNIT_SECOND && (rtv /= 60.0)
                                    rt = rtv
                                end
                            elseif acc == CV_IMS_POSITION_X
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (pos_x = parse(Int, v))
                            elseif acc == CV_IMS_POSITION_Y
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (pos_y = parse(Int, v))
                            elseif acc == CV_IMS_POSITION_Z
                                v = _imzml_attr(l, _IMZML_RE_VALUE); v !== nothing && (pos_z = parse(Int, v))
                            end
                        end
                    end
                end

                # Capture the shared m/z axis (continuous mode) on first encounter
                if is_continuous && !shared_mz_read && !isempty(mz)
                    shared_mz = mz
                    shared_mz_read = true
                end

                # Derive TIC / base peak if the header omitted them (matches DOM reader)
                if tic == 0.0 && !isempty(int_arr)
                    tic = sum(int_arr)
                end
                if basePeakIntensity == 0.0 && !isempty(int_arr)
                    maxidx = argmax(int_arr)
                    basePeakMz = mz[maxidx]
                    basePeakIntensity = int_arr[maxidx]
                end

                metadata = Dict{String,Any}("position_x" => pos_x, "position_y" => pos_y)
                pos_z != 0 && (metadata["position_z"] = pos_z)

                scans[idx] = MSscans(idx, rt, tic, mz, int_arr, msLevel,
                                     basePeakMz, basePeakIntensity, 0.0, polarity,
                                     "", 0.0, 0, spectrumType, -1.0, 0.0, :none, metadata)
                prog === nothing || next!(prog)
            end

            idx == scanCount || throw(ImzmlStreamFallback("found $idx spectra, expected $scanCount"))
        end

        return MSrun(scans)
    finally
        close(io)
    end
end
