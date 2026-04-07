"""
Interface to the mzxml file format
"""



"""
    info_mzxml(filename::String, info::Vector{String}, verbose::Bool=false)
Returns the information content of an mzXML file into a string. Verbosity is controlled by the verbose Boolean variable set by default to false.

"""
function info_mzxml(filename::String, info::Vector{String}, verbose::Bool=false)
    filter = ""
    xdoc = parse_file(filename)
    xroot = root(xdoc)       
    if name(xroot) != "mzXML"
        return ErrorException("Not an mzXML file.")
    end
    
    msRun = find_element(xroot, "msRun")
    scanCount = attribute(msRun, "scanCount")

    for c in child_elements(msRun)     
        if verbose == true
            if name(c) == "parentFile"
                filter = ""
                filter *= name(c) * ": " * attribute(c, "fileName")
                push!(info, filter)
                filter = ""
            elseif name(c) == "msInstrument"       
                for a in child_elements(c)              
                    if has_attribute(a, "category")
                        filter *= attribute(a, "category") * ": " * attribute(a, "value")
                        push!(info, filter)
                        filter = ""
                    end
                    if name(a) == "operator"
                        filter *= name(a) * ": " * attribute(a, "first") * " " * attribute(a, "last")
                        push!(info, filter)
                        filter = ""
                    end 
                    if name(a) == "software"
                        filter *= name(a) * ": " * attribute(a, "name") * ", " * attribute(a, "version")
                        push!(info, filter)
                        filter = ""
                    end
                end  
            elseif name(c) == "dataProcessing"
                for a in child_elements(c)
                    filter *= name(c) * ": " * attribute(a, "type") * ", " * attribute(a, "name") * " " * attribute(a, "version")
                    push!(info, filter)
                    filter = ""
                end        
            end
        end
        while name(c) == "scan"
            filter = scanCount * " scans"
            if !(filter in info)!
                filter !="" ? push!(info, filter) : nothing
            end
            filter = ""

            msLevel = attribute(c, "msLevel")
            filter *=  "MS" * msLevel
            if has_attribute(c, "polarity")
                polarity = attribute(c, "polarity")
                filter *= polarity
            end
            if find_element(c, "precursorMz") != nothing
                precursorMz = find_element(c, "precursorMz")
                filter *= " " * content(precursorMz) * " "

                if has_attribute(precursorMz, "activationMethod")
                    activationMethod = attribute(precursorMz, "activationMethod")
                    filter *= " " * activationMethod
                end
            end
            if has_attribute(c, "collisionEnergy")
                collisionEnergy = attribute(c, "collisionEnergy")
                filter *= "(CE=" * collisionEnergy * ")"
            end
            
            if !(filter in info)!
                filter !="" ? push!(info, filter) : nothing
            end
            filter = ""
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end   

    free(xdoc) ;
    info
end

"""
    load_mzxml_all(filename::String)
Load an entire an mzxml file. 
"""
function load_mzxml_all(filename::String)
    xdoc = parse_file(filename)
    xroot = root(xdoc)
    if name(xroot) != "mzXML"
        return ErrorException("Not an mzXML file.")
    end
    msRun = find_element(xroot, "msRun")
    scanCount = attribute(msRun, "scanCount")
    
    scans = Vector{MSscan}(undef, parse(Int, scanCount))
    index = 1              
    for c1 in child_elements(msRun)
        while name(c1) == "scan"
            scans[index] = load_mzxml_spectrum(c1)
            index += 1
            c1 = find_element(c1,"scan")
            if c1 == nothing
                break
            end
        end
    end
    free(xdoc)   
    return scans
end

"""
    load_mzxml(filename::String, index::Int)
Load from an mzxml file the scan num that matches the input index.
"""
function load_mzxml(filename::String, index::Int)
    xdoc = parse_file(filename)
    xroot = root(xdoc)
    if name(xroot) != "mzXML"
        return ErrorException("Not an mzXML file.")
    end
    msRun = find_element(xroot, "msRun")
    scanCount = attribute(msRun, "scanCount")
    
    for c in child_elements(msRun)
        while name(c) == "scan"
            if parse(Int,attribute(c, "num")) == index
                scan = load_mzxml_spectrum(c)
                free(xdoc)
                return scan
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end   
   
end


"""
    load_mzxml_spectrum(c::XMLElement)
From an XMLElement, returns the data into an MSscan. 
"""
function load_mzxml_spectrum(c::XMLElement)
    num = attribute(c, "num")
    msLevel = attribute(c, "msLevel")
    polarity = ""
    activationMethod = ""
    
    if has_attribute(c, "polarity")
        polarity = attribute(c, "polarity")
    end
    
    if find_element(c, "precursorMz") != nothing
        precursorMz = find_element(c, "precursorMz")
        precursor = content(precursorMz)
        if attribute(precursorMz, "activationMethod") != nothing
            activationMethod = attribute(precursorMz, "activationMethod")
        end
        
    else
        precursor = "0"
        activationMethod = ""
    end

    if has_attribute(c, "basePeakMz",)
        basePeakMz = attribute(c, "basePeakMz")
        basePeakIntensity = attribute(c, "basePeakIntensity")
    else
        basePeakMz = 0.
        basePeakIntensity = 0.
    end
    
    if has_attribute(c, "collisionEnergy")
        collisionEnergy = attribute(c, "collisionEnergy")
    else
        collisionEnergy = "0"
    end
    
    if has_attribute(c, "totIonCurrent")
        totIonCurrent = attribute(c, "totIonCurrent")
        retentionTime = attribute(c, "retentionTime")
    end
    
    peaks = find_element(c, "peaks")
    pairOrder = attribute(peaks, "pairOrder")
    if pairOrder == nothing
        pairOrder = attribute(peaks, "contentType")
    end
    if pairOrder != "m/z-int"
        return MSscan(0 , 0.0, 0.0, [], [], 0, 0.0, 0.0 , 0.0, "", "", 0.0 )
    end

    
    compression = attribute(peaks, "compressionType")
    if compression == "zlib"
        len = attribute(peaks, "compressedLen")       
        data = Libz.inflate( decode( Base64, content(peaks) ) )
    else   
        data = decode(Base64, content(peaks))
    end
    
    precision = attribute(peaks, "precision")
    if precision == "32"
        A = reinterpret(Float32, data)
    elseif precision == "64"
        A = reinterpret(Float64, data)
    end
    
    byteOrder = attribute(peaks, "byteOrder")
    if byteOrder == "network"
        A = ntoh.(A)
    else
        return MSscan(0 , 0.0, 0.0, [], [], 0, 0.0, 0.0 , 0.0, "", "", 0.0 )
    end
    int = convert(Array{Float64,1}, A[2:2:end])
    mz  = convert(Array{Float64,1}, A[1:2:end])
    
    return MSscan(parse(Int,num) , parse(Float64, retentionTime[3:end-1]), parse(Float64,totIonCurrent), mz, int, parse(Int, msLevel), parse(Float64, basePeakMz), parse(Float64, basePeakIntensity), parse(Float64, precursor), polarity, activationMethod, parse(Float64, collisionEnergy) )
    
end

"""
    retention_time(msRun::XMLElement)
From an XMLElement returns the retention time.
"""
function retention_time(msRun::XMLElement)
    rt  = Vector{Float64}(undef,0)
    for c1 in child_elements(msRun)
        while name(c1) == "scan"
            if has_attribute(c1, "totIonCurrent")
                retentionTime = attribute(c1, "retentionTime")
                push!(rt, parse(Float64,retentionTime[3:end-1]))
            end           
            c1 = find_element(c1,"scan")
            if c1 == nothing
                break
            end
        end
    end    
    return rt
end


"""
    filter(msRun::XMLElement, argument::Precursor{<:Real})
Search for scans matching the argument precursor mz and returns a list of index
"""
function filter(msRun::XMLElement, argument::Precursor{<:Real})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).precursor == argument.arg
                push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end


"""
    filter(msRun::XMLElement, argument::Level{<:Int})
Search for scans matching the argument level and returns a list of index

"""
function filter(msRun::XMLElement, argument::Level{<:Int})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).level == argument.arg
                    push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Level{<:AbstractVector}
Search for scans matching the argument levels and returns a list of index

"""
function filter(msRun::XMLElement, argument::Level{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).level == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Precursor{<:AbstractVector})
Search for scans matching the argument precusors mz and returns a list of index

"""
function filter(msRun::XMLElement, argument::Precursor{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).precursor == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end


"""
    filter(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
Search for scans matching the argument activation energies and returns a list of index
"""
function filter(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).collisionEnergy == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Activation_Energy{<:Real})
Search for scans matching the argument activation energy and returns a list of index
"""
function filter(msRun::XMLElement, argument::Activation_Energy{<:Real})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).collisionEnergy == argument.arg
                    push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
Search for scans matching the argument activation methods and returns a list of index
"""
function filter(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).activationMethod == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Activation_Method{<:String})
Search for scans matching the argument activation method and returns a list of index
"""
function filter(msRun::XMLElement, argument::Activation_Method{<:String})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).activationMethod == argument.arg
                    push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Polarity{<:AbstractVector})
Search for scans matching the argument polarities and returns a list of index
"""
function filter(msRun::XMLElement, argument::Polarity{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).polarity == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Polarity{<:String})
Search for scans matching the argument polarity and returns a list of index
"""
function filter(msRun::XMLElement, argument::Polarity{<:String})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).polarity == argument.arg
                    push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end


"""
    filter(msRun::XMLElement, argument::Scan{<:AbstractVector}
Search for scans matching the argument scan nums and returns a list of index
"""
function filter(msRun::XMLElement, argument::Scan{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for c in child_elements(msRun)
            while name(c) == "scan"
                if load_mzxml_spectrum(c).num == i
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::Scan{<:Int})
Search for scans matching the argument scan num and returns a list of index
"""
function filter(msRun::XMLElement, argument::Scan{<:Int})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).num == argument.arg
                push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::RT{<:Real})
Search for scans matching the argument retention time and returns a list of index
"""
function filter(msRun::XMLElement, argument::RT{<:Real})
    subindex = Set{Int}()
    rt = retention_time(msRun)
    index = num2pnt(rt, argument.arg)

    for c in child_elements(msRun)
        while name(c) == "scan"
            if load_mzxml_spectrum(c).num == index
                    push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::RT{<:AbstractVector}
Search for scans matching the argument retention times within the range and returns a list of index
"""
function filter(msRun::XMLElement, argument::RT{<:AbstractVector})
    subindex = Set{Int}()
    rt = retention_time(msRun)
    for i in argument.arg
        index_low  = num2pnt(rt, argument.arg[1])
        index_high = num2pnt(rt, argument.arg[2])
        for c in child_elements(msRun)
            while name(c) == "scan"
                if index_low <= load_mzxml_spectrum(c).num <= index_high
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
Search for scans matching the argument retention times within the ranges and returns a list of index
"""
function filter(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
    subindex = Set{Int}()
    rt = retention_time(msRun)
    for el in argument.arg
        index_low  = num2pnt(rt, el[1])
        index_high = num2pnt(rt, el[2])
        for c in child_elements(msRun)
            while name(c) == "scan"
                if index_low <= load_mzxml_spectrum(c).num <= index_high
                    push!(subindex, load_mzxml_spectrum(c).num)
                end
                c = find_element(c,"scan")
                if c == nothing
                    break
                end
            end
        end
    end
    return subindex
end

"""
    filter(msRun::XMLElement, argument::IC{<:AbstractVector})
Search for scans matching for which the total ion current is within the input range returns a list of index
"""
function filter(msRun::XMLElement, argument::IC{<:AbstractVector})
    subindex = Set{Int}()
    for c in child_elements(msRun)
        while name(c) == "scan"
            if argument.arg[1] <= load_mzxml_spectrum(c).tic <= argument.arg[2]
                push!(subindex, load_mzxml_spectrum(c).num)
            end
            c = find_element(c,"scan")
            if c == nothing
                break
            end
        end
    end
    return subindex
end


"""
    extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
Returns the extracted chromatogram for input file according to the selected method and for set of scan num as input
"""
function extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    if method isa BasePeak
        for i = 1:length(indices)
            push!(xrt, load_mzxml(filename, indices[i]).rt)
            push!(xic, load_mzxml(filename, indices[i]).basePeakIntensity)
        end
    elseif method isa ∆MZ
        mz1 = convert(Float64, method.arg[1] - method.arg[2] )  # mz - ∆mz
        if(mz1 < 0.0)
            return ErrorException("Bad mz ± ∆mz values.")
        end
        mz2 = convert(Float64, method.arg[1] + method.arg[2] ) # mz + ∆mz
        for i = 1:length(indices)
            value = add_ion_current(load_mzxml(filename, indices[i]).mz, load_mzxml(filename, indices[i]).int, mz1, mz2)
            push!(xrt, load_mzxml(filename, indices[i]).rt)
            push!(xic, value)
        end
        
    elseif method isa MZ
        mz1 = convert(Float64,method.arg[1])
        mz2 = convert(Float64,method.arg[2])
        for i = 1:length(indices)
            value = add_ion_current(load_mzxml(filename, indices[i]).mz, load_mzxml(filename, indices[i]).int, mz1, mz2)
            push!(xrt, load_mzxml(filename, indices[i]).rt)
            push!(xic, value)
        end
    else
        for i = 1:length(indices)
            push!(xrt, load_mzxml(filename, indices[i]).rt)
            push!(xic, load_mzxml(filename, indices[i]).tic)
        end

    end
    return Chromatogram(xrt, xic, maximum(xic))
end

"""
    composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
Returns the average MSscans for input filename and according to the input scan num. Calculation of variance is controlled by the stats Boolean variable.
"""
function composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
    if stats == false
        result = load_mzxml(filename, indices[1])
        for i = 2:length(indices)
            result += load_mzxml(filename, indices[i])
        end
        return result / length(indices)
    elseif stats == true
        result = load_mzxml(filename, indices[1])
        for i = 2:length(indices)
            result = avg(result, load_mzxml(filename, indices[i]))
        end
        return standard_deviation(result, length(indices))
    end

end

