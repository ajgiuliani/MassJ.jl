"""
Interface to the MSP (NIST Mass Spectral Library) file format.
MSP is a text-based format used by NIST, MoNA, MassBank, and GNPS for spectral libraries.
Each entry contains metadata fields followed by a peak list.
"""


"""
    load_msp_all(filename::String)
Load all spectra from an MSP file. Returns a `Vector{MSscan}`.
Each entry (Name: ... Num Peaks: ... peak data ... blank line) becomes one MSscan.
"""
function load_msp_all(filename::String)
    scans = Vector{MSscan}(undef, 0)
    scan_index = 0

    open(filename, "r") do io
        params = Dict{String,String}()
        mz_vals = Float64[]
        int_vals = Float64[]
        in_peaks = false
        peaks_remaining = 0

        for raw_line in eachline(io)
            line = strip(raw_line)

            # Blank line = entry separator
            if isempty(line)
                if !isempty(params)
                    scan_index += 1
                    push!(scans, build_msp_scan(params, mz_vals, int_vals, scan_index))
                    empty!(params)
                    empty!(mz_vals)
                    empty!(int_vals)
                    in_peaks = false
                    peaks_remaining = 0
                end
                continue
            end

            # Skip comment lines
            if startswith(line, "#") || startswith(line, ";")
                continue
            end

            if in_peaks && peaks_remaining > 0
                # Parse peak data: handle "mz int; mz int;" and "mz int" formats
                pairs = split(line, ";")
                for pair in pairs
                    pair = strip(pair)
                    isempty(pair) && continue
                    # Split on whitespace, handle optional annotation column
                    parts = split(pair)
                    if length(parts) >= 2
                        mz_val = tryparse(Float64, parts[1])
                        int_val = tryparse(Float64, parts[2])
                        if mz_val !== nothing && int_val !== nothing
                            push!(mz_vals, mz_val)
                            push!(int_vals, int_val)
                            peaks_remaining -= 1
                        end
                    end
                end
            else
                # Metadata line: KEY: value
                colonpos = findfirst(':', line)
                if colonpos !== nothing
                    key = uppercase(strip(line[1:colonpos-1]))
                    value = strip(line[colonpos+1:end])

                    if key == "NUM PEAKS"
                        params[key] = value
                        np = tryparse(Int, value)
                        if np !== nothing && np > 0
                            in_peaks = true
                            peaks_remaining = np
                        end
                    else
                        # Accumulate repeated fields (e.g., Synon)
                        if haskey(params, key) && key in ("SYNON", "COMMENT", "COMMENTS")
                            params[key] *= "; " * value
                        else
                            params[key] = value
                        end
                    end
                else
                    # Could be a peak line without prior "Num Peaks" (non-standard)
                    # or continuation of peaks
                    if in_peaks
                        parts = split(line)
                        if length(parts) >= 2
                            mz_val = tryparse(Float64, parts[1])
                            int_val = tryparse(Float64, parts[2])
                            if mz_val !== nothing && int_val !== nothing
                                push!(mz_vals, mz_val)
                                push!(int_vals, int_val)
                                peaks_remaining -= 1
                            end
                        end
                    end
                end
            end
        end

        # Handle last entry (no trailing blank line)
        if !isempty(params)
            scan_index += 1
            push!(scans, build_msp_scan(params, mz_vals, int_vals, scan_index))
        end
    end

    return scans
end


"""
    build_msp_scan(params::Dict{String,String}, mz::Vector{Float64}, int::Vector{Float64}, index::Int)
Build an MSscan from parsed MSP entry parameters and peak data.
"""
function build_msp_scan(params::Dict{String,String}, mz::Vector{Float64},
                        int::Vector{Float64}, index::Int)
    mz_copy = copy(mz)
    int_copy = copy(int)

    num = index

    # Retention time
    rt = 0.0
    for key in ("RT", "RETENTION_TIME", "RETENTIONTIME")
        if haskey(params, key)
            r = tryparse(Float64, split(params[key])[1])
            if r !== nothing
                rt = r
                break
            end
        end
    end

    # TIC
    tic = isempty(int_copy) ? 0.0 : sum(int_copy)

    # MS level
    msLevel = 2  # MSP is typically MS/MS
    for key in ("SPECTRUM_TYPE", "MSLEVEL", "MS_LEVEL")
        if haskey(params, key)
            val = uppercase(strip(params[key]))
            if val == "MS1" || val == "1"
                msLevel = 1
            elseif startswith(val, "MS") && length(val) > 2
                ml = tryparse(Int, val[3:end])
                if ml !== nothing
                    msLevel = ml
                end
            else
                ml = tryparse(Int, val)
                if ml !== nothing
                    msLevel = ml
                end
            end
            break
        end
    end

    # Base peak
    basePeakMz = 0.0
    basePeakIntensity = 0.0
    if !isempty(int_copy)
        maxidx = argmax(int_copy)
        basePeakMz = mz_copy[maxidx]
        basePeakIntensity = int_copy[maxidx]
    end

    # Precursor m/z
    precursorMz = 0.0
    for key in ("PRECURSOR_MZ", "PRECURSORMZ", "PEPMASS")
        if haskey(params, key)
            p = tryparse(Float64, split(params[key])[1])
            if p !== nothing
                precursorMz = p
                break
            end
        end
    end

    # Polarity / ion mode
    polarity = ""
    for key in ("ION_MODE", "IONMODE", "ION MODE")
        if haskey(params, key)
            val = uppercase(strip(params[key]))
            if startswith(val, "POS") || val == "P"
                polarity = "+"
            elseif startswith(val, "NEG") || val == "N"
                polarity = "-"
            end
            break
        end
    end

    # Collision energy
    collisionEnergy = 0.0
    for key in ("COLLISION_ENERGY", "COLLISIONENERGY", "CE")
        if haskey(params, key)
            # Extract numeric part (e.g., "35 eV" -> 35.0)
            ce = tryparse(Float64, split(params[key])[1])
            if ce !== nothing
                collisionEnergy = ce
                break
            end
        end
    end

    # Charge state from precursor type (e.g., "[M+H]+", "[M+2H]2+")
    chargeState = 0

    # Activation method
    activationMethod = ""

    # Metadata
    metadata = Dict{String,Any}()
    if haskey(params, "NAME")
        metadata["name"] = params["NAME"]
    end
    if haskey(params, "FORMULA")
        metadata["formula"] = params["FORMULA"]
    end
    if haskey(params, "MW")
        metadata["mw"] = params["MW"]
    end
    if haskey(params, "INCHIKEY")
        metadata["inchikey"] = params["INCHIKEY"]
    end
    if haskey(params, "PRECURSOR_TYPE") || haskey(params, "PRECURSORTYPE")
        ptype = get(params, "PRECURSOR_TYPE", get(params, "PRECURSORTYPE", ""))
        metadata["precursor_type"] = ptype
    end
    if haskey(params, "COMMENT") || haskey(params, "COMMENTS")
        metadata["comments"] = get(params, "COMMENTS", get(params, "COMMENT", ""))
    end
    if haskey(params, "SYNON")
        metadata["synonyms"] = params["SYNON"]
    end
    if haskey(params, "CAS#")
        metadata["cas"] = params["CAS#"]
    end
    if haskey(params, "DB#")
        metadata["db_id"] = params["DB#"]
    end
    if haskey(params, "NIST#")
        metadata["nist_id"] = params["NIST#"]
    end
    if haskey(params, "SPLASH")
        metadata["splash"] = params["SPLASH"]
    end
    if haskey(params, "INSTRUMENT")
        metadata["instrument"] = params["INSTRUMENT"]
    end
    if haskey(params, "INSTRUMENT_TYPE")
        metadata["instrument_type"] = params["INSTRUMENT_TYPE"]
    end

    return MSscan(num, rt, tic, mz_copy, int_copy, msLevel,
                  basePeakMz, basePeakIntensity, precursorMz, polarity,
                  activationMethod, collisionEnergy,
                  chargeState, :centroid, -1.0, 0.0, :none, metadata)
end


"""
    info_msp(filename::String, info::Vector{String}, verbose::Bool=false)
Returns summary information about an MSP file.
"""
function info_msp(filename::String, info::Vector{String}, verbose::Bool=false)
    scan_count = 0
    seen = Set{String}()

    open(filename, "r") do io
        params = Dict{String,String}()

        for raw_line in eachline(io)
            line = strip(raw_line)

            if isempty(line)
                if !isempty(params)
                    scan_count += 1

                    # Build description
                    msLevel = "2"
                    for key in ("SPECTRUM_TYPE", "MSLEVEL", "MS_LEVEL")
                        if haskey(params, key)
                            val = uppercase(strip(params[key]))
                            if startswith(val, "MS")
                                msLevel = val[3:end]
                            else
                                msLevel = val
                            end
                            break
                        end
                    end
                    desc = "MS" * msLevel

                    for key in ("ION_MODE", "IONMODE")
                        if haskey(params, key)
                            val = uppercase(strip(params[key]))
                            if startswith(val, "POS")
                                desc *= "+"
                            elseif startswith(val, "NEG")
                                desc *= "-"
                            end
                            break
                        end
                    end

                    for key in ("PRECURSOR_MZ", "PRECURSORMZ")
                        if haskey(params, key)
                            desc *= " " * split(params[key])[1]
                            break
                        end
                    end

                    if !(desc in seen)
                        push!(seen, desc)
                        push!(info, desc)
                    end

                    if verbose && haskey(params, "NAME")
                        entry = "Name: " * params["NAME"]
                        if !(entry in seen)
                            push!(seen, entry)
                            push!(info, entry)
                        end
                    end

                    empty!(params)
                end
                continue
            end

            colonpos = findfirst(':', line)
            if colonpos !== nothing
                key = uppercase(strip(line[1:colonpos-1]))
                value = strip(line[colonpos+1:end])
                params[key] = value
            end
        end

        # Handle last entry
        if !isempty(params)
            scan_count += 1
            msLevel = "2"
            for key in ("SPECTRUM_TYPE", "MSLEVEL", "MS_LEVEL")
                if haskey(params, key)
                    val = uppercase(strip(params[key]))
                    if startswith(val, "MS")
                        msLevel = val[3:end]
                    else
                        msLevel = val
                    end
                    break
                end
            end
            desc = "MS" * msLevel

            for key in ("ION_MODE", "IONMODE")
                if haskey(params, key)
                    val = uppercase(strip(params[key]))
                    if startswith(val, "POS")
                        desc *= "+"
                    elseif startswith(val, "NEG")
                        desc *= "-"
                    end
                    break
                end
            end

            for key in ("PRECURSOR_MZ", "PRECURSORMZ")
                if haskey(params, key)
                    desc *= " " * split(params[key])[1]
                    break
                end
            end

            if !(desc in seen)
                push!(seen, desc)
                push!(info, desc)
            end

            empty!(params)
        end
    end

    pushfirst!(info, string(scan_count) * " scans")
    return info
end
