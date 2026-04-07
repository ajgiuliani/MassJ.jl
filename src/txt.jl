"""
Interface to txt file format
"""

"""
    info_txt!(filename::String, info::Vector{String})
Returns the information content of a txt file into a string.

"""
function info_txt!(filename::String, info::Vector{String})
    counter = 1
    filter = ""
"""
    for (index,lines) in enumerate(readlines(filename))        
        rt = 0.0
        polarity = ""
        ionizationMethod=""
        activationMethod = ""
        collisionEnergy = 0.0
        msLevel=0
        precursorMz=0.0
        points=0
        for (index, el) in enumerate(split(lines, ","))
            if index == 1
                #rt = parse.(Float64, el)
            elseif index == 2
                polarity = el
                filter *= polarity
            elseif index == 3
                ionizationMethod = el
            elseif index == 4
                msLevel = parse.(Int, split(el, "ms")[2])
                filter *= "MS" * split(el, "ms")[2]
            elseif index == 5
                precursorMz = parse.(Float64, el)
                filter *= " " * el * " "
            elseif index == 8
                points = parse.(Int, split(el, " ")[1])
            elseif index >= 9
            end
        end
        if !(filter in info)
            filter != "" ? push!(info, filter) : nothing
        end
        filter = ""
        counter +=1
    end
    filter = string(counter) * " scans"
    pushfirst!(info, filter)
"""
    return info
end



"""
    load_txt_all(filename::String)
Load the content of an ascii file into a Vector of MSscan.
# Examples
```julia-repl
julia> scan = msJ.load("test1.txt")
MSscan(10, 11.0422, 5.37656e6, [140.083  …  2000.0], [91.044 … 29.5873])
```
"""
function load_txt_all(filename::String)
    return  load_txt_spectrum(filename)
end


function load_txt_spectrum(filename::String)

    scan = readdlm(filename, Float64)
    
    num = 1
    rt = 0.0
    tic = sum(scan[:,2])
    
    msLevel = 0
    basePeakIntensity = maximum(scan[:,2])
    basePeakMz = scan[:,1][ num2pnt(scan[:,2], basePeakIntensity)]
    precursorMz = 0.0
    polarity = ""
    activationMethod = ""
    collisionEnergy = 0.0
    
    return  MSscan(num, rt, tic, scan[:,1], scan[:,2], msLevel, basePeakMz, basePeakIntensity, precursorMz, polarity, activationMethod, collisionEnergy)
end


