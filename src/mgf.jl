"""
Interface to the MGF (Mascot Generic Format) file format.
MGF is a text-based peak list format that stores centroided MS/MS data.
"""


"""
    load_mgf_all(filename::String)
Load all spectra from an MGF file. Returns a Vector{MSscan}.
Each BEGIN IONS...END IONS block becomes one MSscan.
"""
function load_mgf_all(filename::String)
    scans = Vector{MSscan}(undef, 0)
    scan_index = 0

    open(filename, "r") do io
        in_ions = false
        mz_vals = Float64[]
        int_vals = Float64[]
        params = Dict{String,String}()

        for raw_line in eachline(io)
            line = strip(raw_line)

            if isempty(line) || startswith(line, "#") || startswith(line, ";")
                continue
            end

            uline = uppercase(line)

            if uline == "BEGIN IONS"
                in_ions = true
                empty!(mz_vals)
                empty!(int_vals)
                empty!(params)
                continue
            end

            if uline == "END IONS"
                in_ions = false
                scan_index += 1
                scan = build_mgf_scan(params, mz_vals, int_vals, scan_index)
                push!(scans, scan)
                continue
            end

            if !in_ions
                continue  # skip global parameters
            end

            # Inside an ion block: either KEY=value or peak data
            eqpos = findfirst('=', line)
            if eqpos !== nothing
                key = uppercase(line[1:eqpos-1])
                value = line[eqpos+1:end]
                params[key] = value
            else
                # Peak data line: m/z intensity [charge]
                parts = split(line)
                if length(parts) >= 2
                    mz_val = tryparse(Float64, parts[1])
                    int_val = tryparse(Float64, parts[2])
                    if mz_val !== nothing && int_val !== nothing
                        push!(mz_vals, mz_val)
                        push!(int_vals, int_val)
                    end
                end
            end
        end
    end

    return scans
end


"""
    build_mgf_scan(params::Dict{String,String}, mz::Vector{Float64}, int::Vector{Float64}, index::Int)
Build an MSscan from parsed MGF ion block parameters and peak data.
"""
function build_mgf_scan(params::Dict{String,String}, mz::Vector{Float64},
                        int::Vector{Float64}, index::Int)
    mz_copy = copy(mz)
    int_copy = copy(int)

    # Scan number: use sequential index to match array position
    # (the filtering system requires scan num == array index)
    num = index

    # Retention time (RTINSECONDS -> convert to minutes)
    rt = 0.0
    if haskey(params, "RTINSECONDS")
        rtstr = split(params["RTINSECONDS"], "-")[1]  # handle ranges
        r = tryparse(Float64, rtstr)
        if r !== nothing
            rt = r / 60.0  # convert to minutes
        end
    end

    # TIC (sum of intensities)
    tic = isempty(int_copy) ? 0.0 : sum(int_copy)

    # MS level (MGF is typically MS2, but check for MSLEVEL parameter)
    msLevel = 2
    if haskey(params, "MSLEVEL")
        ml = tryparse(Int, params["MSLEVEL"])
        if ml !== nothing
            msLevel = ml
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

    # Precursor m/z and intensity
    precursorMz = 0.0
    if haskey(params, "PEPMASS")
        parts = split(params["PEPMASS"])
        p = tryparse(Float64, parts[1])
        if p !== nothing
            precursorMz = p
        end
    end

    # Charge state
    chargeState = 0
    if haskey(params, "CHARGE")
        chargeStr = strip(params["CHARGE"])
        # Handle "2+", "3-", or "2"
        # Take first charge if multiple are given (e.g., "2+ and 3+")
        firstCharge = split(chargeStr, r"\s+and\s+|,")[1]
        firstCharge = strip(firstCharge)
        # Remove sign suffix
        if endswith(firstCharge, "+") || endswith(firstCharge, "-")
            cs = tryparse(Int, firstCharge[1:end-1])
        else
            cs = tryparse(Int, firstCharge)
        end
        if cs !== nothing
            chargeState = cs
        end
    end

    # Polarity (from charge sign)
    polarity = ""
    if haskey(params, "CHARGE")
        chargeStr = strip(params["CHARGE"])
        if contains(chargeStr, "-")
            polarity = "-"
        elseif contains(chargeStr, "+")
            polarity = "+"
        end
    end

    # Activation method
    activationMethod = ""
    collisionEnergy = 0.0

    # Store extra MGF parameters in metadata
    metadata = Dict{String,Any}()
    if haskey(params, "TITLE")
        metadata["title"] = params["TITLE"]
    end
    if haskey(params, "SCANS")
        metadata["scans"] = params["SCANS"]
    end

    return MSscan(num, rt, tic, mz_copy, int_copy, msLevel,
                  basePeakMz, basePeakIntensity, precursorMz, polarity,
                  activationMethod, collisionEnergy,
                  chargeState, :centroid, -1.0, 0.0, :none, metadata)
end


"""
    info_mgf(filename::String, info::Vector{String}, verbose::Bool=false)
Returns summary information about an MGF file.
"""
function info_mgf(filename::String, info::Vector{String}, verbose::Bool=false)
    scan_count = 0
    seen = Set{String}()

    open(filename, "r") do io
        in_ions = false
        params = Dict{String,String}()

        for raw_line in eachline(io)
            line = strip(raw_line)
            uline = uppercase(line)

            if uline == "BEGIN IONS"
                in_ions = true
                empty!(params)
                continue
            end

            if uline == "END IONS"
                in_ions = false
                scan_count += 1

                # Build description
                msLevel = "2"
                if haskey(params, "MSLEVEL")
                    msLevel = params["MSLEVEL"]
                end
                desc = "MS" * msLevel

                if haskey(params, "CHARGE")
                    chargeStr = strip(params["CHARGE"])
                    if contains(chargeStr, "-")
                        desc *= "-"
                    elseif contains(chargeStr, "+")
                        desc *= "+"
                    end
                end

                if haskey(params, "PEPMASS")
                    pmz = split(params["PEPMASS"])[1]
                    desc *= " " * pmz
                end

                if !(desc in seen)
                    push!(seen, desc)
                    push!(info, desc)
                end
                continue
            end

            if in_ions
                eqpos = findfirst('=', line)
                if eqpos !== nothing
                    key = uppercase(line[1:eqpos-1])
                    value = line[eqpos+1:end]
                    params[key] = value
                end
            end
        end
    end

    pushfirst!(info, string(scan_count) * " scans")
    return info
end
