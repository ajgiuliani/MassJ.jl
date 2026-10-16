"""
Module for extracting subsets from a Vector{MSscans} according to specific conditions
"""

# User Interface.
# ---------------

export extract



"""
    extract(filename::String, arguments::FilterType...)
Returns a `Vector{MSscans}` containing the scans that match the given [`FilterType`](@ref) conditions. Without arguments, returns all scans.

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

# Examples
```julia-repl
julia> sub_set = extract("test.mzXML", MassJ.Level(2))
2-element Vector{MassJ.MSscans}:
 MSscans(num=2, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=0.7307 min, tic=9727.2  precursor=1255.5)
 MSscans(num=5, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=4.3442 min, tic=12203.5  precursor=1255.5)

julia> sub_set = extract("test.mzML", MassJ.Level(1))
1-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS1+, 5 pts m/z=[100.0, 500.0], rt=0.5 min, tic=19000.0)

julia> sub_set = extract("test.msp", MassJ.Polarity("+"))
2-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS2+@20.0eV, 5 pts m/z=[42.034, 195.088], rt=0.0 min, tic=178600.0  precursor=195.0877)
 MSscans(num=3, MS1+, 4 pts m/z=[180.063, 183.074], rt=0.0 min, tic=8955.0)
```
"""
function extract(filename::String, arguments::FilterType...)
    extension = split(filename, ".")[end]
    ext = Unicode.normalize(extension, casefold=true)

    if ext == "mzxml"
        index = Set{Int}()
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            error("Not an mzXML file.")
        end
        msRun = find_element(xroot, "msRun")
        scanCount = parse(Int, attribute(msRun, "scanCount"))

        index = Set( i for i in 1:scanCount )

        for el in arguments
            subindex = filter(msRun, el)
            index = intersect(index, subindex)
        end

        free(xdoc)
        indices = sort([ i for i in index])
        if length(indices) >= 1
            return build_subset(filename, indices)
        else
            error("No matching spectra.")
        end

    elseif ext in ("mzml", "mgf", "msp", "imzml")
        scans = load(filename)
        return extract(scans, arguments...)

    else
        error("File format not supported.")
    end
end



"""
    build_subset(filename::String, indices::Vector{Int})
Returns a Vector of `MSscans` from the input file according to the scan num (indices).
"""
function build_subset(filename::String, indices::Vector{Int})
    sub_set = Vector{MSscans}(undef,0)   
    for i = 1:length(indices)
        push!(sub_set, load_mzxml(filename, indices[i]))
    end
    return sub_set
end


"""
    extract(scans::AbstractVector{MSscans}, arguments::FilterType...)
Search for scans matching the argument MS level and returns an array of matching MSscans; throws an error ("No matching spectra found.") when nothing matches.
# Examples
```julia-repl
julia> scans = load("test.mzXML");                         # load mass spectra

julia> sub_set = extract(scans)                            # no conditions → original data
6-element Vector{MassJ.MSscans}:
 MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=5.08195e6)
 MSscans(num=2, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=0.7307 min, tic=9727.2  precursor=1255.5)
 MSscans(num=3, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=2.1379 min, tic=11.3032  precursor=902.33)
 MSscans(num=4, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=3.7578 min, tic=4.8084e6)
 MSscans(num=5, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=4.3442 min, tic=12203.5  precursor=1255.5)
 MSscans(num=6, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=5.7689 min, tic=4.84455  precursor=902.33)

julia> sub_set = extract(scans, MassJ.Level(2))            # extract MS/MS spectra
2-element Vector{MassJ.MSscans}:
 MSscans(num=2, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=0.7307 min, tic=9727.2  precursor=1255.5)
 MSscans(num=5, MS2+ CID@18.0eV, 19860 pts m/z=[345.083, 2000.0], rt=4.3442 min, tic=12203.5  precursor=1255.5)
```

"""
function extract(scans::AbstractVector{MSscans}, arguments::FilterType...)
    pred = compose_predicates(scans, arguments)

    sub_set = Vector{MSscans}(undef, 0)
    for scan in scans
        pred(scan) || continue
        push!(sub_set, scan)
    end

    if !isempty(sub_set)
        return sub_set
    else
        error("No matching spectra found.")
    end
end

"""
    build_subset(scans::AbstractVector{MSscans}, indices::Vector{Int})
Returns a Vector of `MSscans` according to the input scan num.
"""
function build_subset(scans::AbstractVector{MSscans}, indices::Vector{Int})
    sub_set = Vector{MSscans}(undef,0)
    for i = 1:length(indices)
        push!(sub_set, scans[indices[i]])
    end
    return sub_set
end



# `extract`/`chromatogram`/`average`/`retention_time` accept any
# `AbstractVector{MSscans}`, so an `MSrun` (which is one) dispatches to those
# methods directly — no separate adapter methods are needed.
