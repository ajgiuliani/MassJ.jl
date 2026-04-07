"""
Module for importing and exporting data. Dispatch to specific methods according to the file extension.
"""


# User Interface.
# ---------------

export info, load, chromatogram, average



"""
    info(filename::String; verbose::Bool = false)
Reads the content of a mass spectrometry file and returns a `Vector{String}` containing the number of scans and the different scan types described by their MS level, polarity and, for MS/MS data, the precursor m/z followed by the activation method and collision energy. Each entry is unique, providing a summary of the input file. With `verbose = true`, the function also returns instrument metadata (parentFile, manufacturer, model, ionisation, mass analyzer, detector, software, data processing) when available.

Supported file formats: mzXML, mzML, MGF.

# Examples
```julia-repl
julia> info("test.mzXML")
4-element Array{String,1}:
 "51 scans"
 "MS1+"
 "MS2+ 1255.5  CID(CE=18)"
 "MS3+ 902.33  PQD(CE=35)"

julia> info("test.mzML")
4-element Array{String,1}:
 "3 scans"
 "MS1+"
 "MS2+ 400.0  CID(CE=25.0)"
 "MS2+ 500.0  HCD(CE=30.0)"
```
"""
function info(filename::String; verbose::Bool = false)
    info = Vector{String}(undef,0)
    extension = split(filename, ".")[end]

    ext = Unicode.normalize(extension, casefold=true)
    if ext == "mzxml"
        return info_mzxml(filename, info, verbose)
    elseif ext == "mzml"
        return info_mzml(filename, info, verbose)
    elseif ext == "mgf"
        return info_mgf(filename, info, verbose)
    else
        return ErrorException("File format not supported.")
    end
end



"""
    load(filename::String)
Loads the mass spectra from a file and returns a `Vector{MSscan}` where each element contains one scan. The function dispatches to the appropriate reader based on the file extension.

Supported file formats: mzXML, mzML, MGF, TXT.

# Examples
```julia-repl
julia> scans = load("test.mzXML")
6-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, ...)

julia> scans = load("test.mzML")
3-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.5, 19000.0, ...)

julia> scans = load("test.mgf")
3-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.5, 4800.0, ...)
```
"""
function load(filename::String)
    extension = split(filename, ".")[end]

    ext = Unicode.normalize(extension, casefold=true)
    if ext == "mzxml"
        return load_mzxml_all(filename)
    elseif ext == "mzml"
        return load_mzml_all(filename)
    elseif ext == "mgf"
        return load_mgf_all(filename)
    elseif ext == "txt"
        return load_txt_all(filename)
    else
        return ErrorException("File format not supported.")
    end
end



"""
    retention_time(filename::String)
Returns a `Vector{Float64}` with the retention times (in minutes) of each scan in the file.

Supported file formats: mzXML, mzML, MGF.

# Examples
```julia-repl
julia> retention_time("test.mzXML")
6-element Array{Float64,1}:
  0.1384
  0.7307
  2.1379
  ...
```
"""
function retention_time(filename::String)
    extension = split(filename, ".")[end]
    rt  = Vector{Float64}(undef,0)

    ext = Unicode.normalize(extension, casefold=true)
    if ext == "mzxml"
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
        end

        msRun = find_element(xroot, "msRun")
        scanCount = attribute(msRun, "scanCount")

        rt = retention_time(msRun)
        free(xdoc)
    elseif ext == "mzml"
        rt = retention_time_mzml(filename)
    elseif ext == "mgf"
        scans = load_mgf_all(filename)
        rt = [s.rt for s in scans]
    else
        ErrorException("File format not supported.")
    end
    return rt
end


"""
    retention_time(scans::Vector{MSscan})
Returns an array composed of the retention times of the individual mass spectra. 
# Examples
```julia-repl
julia> retention_time("scans")
51-element Array{Float64,1}:
  0.1384
  0.7307
  2.1379
....
```
"""
function retention_time(scans::Vector{MSscan})
    rt  = Vector{Float64}(undef,0)
    for elem in scans
        push!(rt, elem.rt)
    end
    return rt
end



"""
    chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
Returns a [`Chromatogram`](@ref) holding the retention time (rt), the ion current (ic) and the maximum value (maxic) for all the mass spectra within the file. The `method` keyword controls how the ion current is computed: `TIC()` (default) for total ion current, `BasePeak()` for base peak intensity, `∆MZ([mz, Δ])` for extracted ion current around mz ± Δ, or `MZ([mz1, mz2])` for a m/z range.

The data may be filtered using [`FilterType`](@ref) arguments: `Level(N)`, `Precursor(mz)`, `Activation_Method("method")`, `Activation_Energy(ce)`, `Polarity("+")`, `Scan(n)`, `RT([t1,t2])`, `IC([min,max])`.

Supported file formats: mzXML, mzML, MGF.

# Examples
```julia-repl
julia> chrom = chromatogram("test.mzXML")
MSj.Chromatogram([0.1384, ...], [5.08195e6, ...], 5.08195e6)

julia> chrom = chromatogram("test.mzML", method = MSj.BasePeak())
MSj.Chromatogram([0.5, 1.0, 1.5], [8000.0, 2000.0, 1200.0], 8000.0)

julia> chrom = chromatogram("test.mzXML", MSj.Level(1), method = MSj.∆MZ([500,5]))
```
"""
function chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
    extension = split(filename, ".")[end]
    ext = Unicode.normalize(extension, casefold=true)

    if ext == "mzxml"
        xrt = Vector{Float64}(undef,0)
        xic = Vector{Float64}(undef,0)
        index = Set{Int}()
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
        end

        msRun = find_element(xroot, "msRun")
        scanCount = parse(Int, attribute(msRun, "scanCount"))

        index = Set( i for i in 1:scanCount )

        for el in filters
            subindex = filter(msRun, el)
            index = intersect(index, subindex)
        end

        free(xdoc)

        indices = sort([ i for i in index])
        if length(indices) != 0
            return extracted_chromatogram(filename, indices, method)
        else
            ErrorException("No matching spectra.")
        end

    elseif ext == "mzml" || ext == "mgf"
        # Load all scans and delegate to the Vector{MSscan} method
        scans = load(filename)
        return chromatogram(scans, filters...; method=method)

    else
        ErrorException("File format not supported.")
    end
end

"""
    chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
Returns the retention time and the total ion current by default for all the mass spectra within the Array of mass spectrum container MSscan. Alternatively, other options may be supplied such as method = MSj.BasePeak, which returs the base peak intensity, method = MSj.∆MZ([500,5]), which returns the ion current for the range mz = 500 ± 5, or method = MSj.MZ([200,1000]) which return the ion current in the range from m/z 200 to m/z 1000.  The data may be filtered by ms level, precursor mass, activation methods, etc, using the arguments MSj.Level(N), MSj.Precursor(mz), MSj.Activation_Method("method")...
# Examples
```julia-repl
julia> rt, ic = chromatogram("test.mzxml")
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
julia> rt, ic = chromatogram("test.mzxml", method = MSj.BasePeak() )
([0.1384  …  60.4793], [102558.0  …  1.23181])
julia> rt, ic = chromatogram("test.mzxml", method = MSj.∆MZ([500,5]) )
([0.1384  …  60.4793], [46036.6  …  14.2529])
julia> rt, ic = chromatogram("test.mzxml", method = MSj.MZ([200,1000]))
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
```
"""
function chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
    # Ranges of mz value used to compute the tic from
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    index = Set( i for i in 1:length(scans) )

    for el in filters
        subindex = filter(scans, el)
        index = intersect(index, subindex)
    end

    indices = sort([ i for i in index])
    if length(indices) != 0
        return extracted_chromatogram(scans, indices, method)
    else
        ErrorException("No matching spectra.")
    end    
end





"""
    average(filename::String, arguments::FilterType...; stats::Bool=true)
Returns the average mass spectrum as an [`MSscans`](@ref) container, along with the sample standard deviation of the intensities when `stats=true` (default). The data may be filtered using [`FilterType`](@ref) arguments.

Supported file formats: mzXML, mzML, MGF.

# Examples
```julia-repl
julia> spectrum = average("test.mzXML")
MSj.MSscans([1, 2, 3, ...

julia> spectrum = average("test.mzML", MSj.Level(2))
MSj.MSscans([2, 3], ...

julia> spectrum = average("test.mzXML", MSj.Precursor(1255.5), MSj.RT([1, 60]))
MSj.MSscans([2, 5, 8, 11], ...
```
"""
function average(filename::String, arguments::FilterType...; stats::Bool=true)
    extension = split(filename, ".")[end]
    ext = Unicode.normalize(extension, casefold=true)

    if ext == "mzxml"
        index = Set{Int}()
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
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
        if length(indices) >= 2
            return composite_spectra(filename, indices, stats)
        elseif length(indices) == 1
            return load_mzxml(filename, indices[1])
        else
            ErrorException("No matching spectra.")
        end

    elseif ext == "mzml" || ext == "mgf"
        scans = load(filename)
        return average(scans, arguments...; stats=stats)

    else
        ErrorException("File format not supported.")
    end
end


"""
    average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
Returns the average mass spectrum container (MSscans) along with the sample standard deviation of the intensities with stats=true (default) for all the mass spectra within the Array of mass spectrum container MSscan.. The data may be filtered by level, precursor mass, activation methods, etc, using the arguments MSj.Level(N), MSj.Precursor(mz), MSj.Activation_Method("method"), or any combination of these arguments.
# Examples
```julia-repl
julia> spectrum = average("test.mzxml")
MSj.MSscans([1, 2, 3 ....
julia> spectrum = average("test.mzxml", MSj.Level(1) )
MSj.MSscans([1, 4, 7, 10,
julia> spectrum = average("test.mzxml", MSj.Precursor(1255.5) )
MSj.MSscans([2, 5, 8, 11, ...
julia> spectrum = average("test.mzxml", MSj.Activation_Method("PQD") )
MSj.MSscans([3, 6, 9, 12, 15,
julia> spectrum = average("test.mzxml", MSj.Activation_Method("PQD"), MSj.Polarity("+"), MSj.RT([10,20]))
MSj.MSscans([9, 12, 15, 18], ...
```
"""
function average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
    index = Set( i for i in 1:length(scans) )
    
    for el in arguments
        subindex = filter(scans, el)
        index = intersect(index, subindex)
    end

    indices = sort([ i for i in index])
    if length(indices) >= 2
        return composite_spectra(scans, indices, stats)
    elseif length(indices) == 1
        return scans[indices[1] ]
    else
        ErrorException("No matching spectra.")
    end
end

