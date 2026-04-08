"""
Module for importing and exporting data. Dispatch to specific methods according to the file extension.
"""


# User Interface.
# ---------------

export info, load, chromatogram, average



"""
    info(filename::String; verbose::Bool = false)
Reads the content of a mass spectrometry file and returns a `Vector{String}` containing the number of scans and the different scan types described by their MS level, polarity and, for MS/MS data, the precursor m/z followed by the activation method and collision energy. Each entry is unique, providing a summary of the input file. With `verbose = true`, the function also returns instrument metadata (parentFile, manufacturer, model, ionisation, mass analyzer, detector, software, data processing) when available.

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

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

julia> info("library.msp")
4-element Array{String,1}:
 "3 scans"
 "MS2+ 195.0877"
 "MS2- 179.0344"
 "MS1+"
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
    elseif ext == "msp"
        return info_msp(filename, info, verbose)
    elseif ext == "imzml"
        return info_imzml(filename, info, verbose)
    else
        return ErrorException("File format not supported.")
    end
end



"""
    load(filename::String)
Loads the mass spectra from a file and returns a `Vector{MSscan}` where each element contains one scan. The function dispatches to the appropriate reader based on the file extension.

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

# Examples
```julia-repl
julia> scans = load("test.mzXML")
6-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.1384, 5.08195e6, ...)

julia> scans = load("test.mzML")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.5, 19000.0, ...)

julia> scans = load("test.mgf")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.5, 4800.0, ...)

julia> scans = load("library.msp")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.0, 178600.0, ...)

julia> scans = load("sample.imzML")
10000-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.0, 8000.0, ...)
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
    elseif ext == "msp"
        return load_msp_all(filename)
    elseif ext == "imzml"
        return load_imzml_all(filename)
    elseif ext == "txt"
        return load_txt_all(filename)
    else
        return ErrorException("File format not supported.")
    end
end



"""
    retention_time(filename::String)
Returns a `Vector{Float64}` with the retention times (in minutes) of each scan in the file.

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

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
    elseif ext == "mgf" || ext == "msp" || ext == "imzml"
        scans = load(filename)
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

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

# Examples
```julia-repl
julia> chrom = chromatogram("test.mzXML")
MassJ.Chromatogram([0.1384, ...], [5.08195e6, ...], 5.08195e6)

julia> chrom = chromatogram("test.mzML", method = MassJ.BasePeak())
MassJ.Chromatogram([0.5, 1.0, 1.5], [8000.0, 2000.0, 1200.0], 8000.0)

julia> chrom = chromatogram("test.mzXML", MassJ.Level(1), method = MassJ.∆MZ([500,5]))
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

    elseif ext in ("mzml", "mgf", "msp", "imzml")
        # Load all scans and delegate to the Vector{MSscan} method
        scans = load(filename)
        return chromatogram(scans, filters...; method=method)

    else
        ErrorException("File format not supported.")
    end
end

"""
    chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
Returns the retention time and the total ion current by default for all the mass spectra within the Array of mass spectrum container MSscan. Alternatively, other options may be supplied such as method = MassJ.BasePeak, which returs the base peak intensity, method = MassJ.∆MZ([500,5]), which returns the ion current for the range mz = 500 ± 5, or method = MassJ.MZ([200,1000]) which return the ion current in the range from m/z 200 to m/z 1000.  The data may be filtered by ms level, precursor mass, activation methods, etc, using the arguments MassJ.Level(N), MassJ.Precursor(mz), MassJ.Activation_Method("method")...
# Examples
```julia-repl
julia> rt, ic = chromatogram("test.mzxml")
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
julia> rt, ic = chromatogram("test.mzxml", method = MassJ.BasePeak() )
([0.1384  …  60.4793], [102558.0  …  1.23181])
julia> rt, ic = chromatogram("test.mzxml", method = MassJ.∆MZ([500,5]) )
([0.1384  …  60.4793], [46036.6  …  14.2529])
julia> rt, ic = chromatogram("test.mzxml", method = MassJ.MZ([200,1000]))
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
```
"""
function chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
    pred = compose_predicates(scans, filters)

    xrt = Vector{Float64}(undef, 0)
    xic = Vector{Float64}(undef, 0)

    for scan in scans
        pred(scan) || continue
        push!(xrt, scan.rt)
        if method isa BasePeak
            push!(xic, scan.basePeakIntensity)
        elseif method isa ∆MZ
            mz1 = convert(Float64, method.arg[1] - method.arg[2])
            if mz1 < 0.0
                return ErrorException("Bad mz ± ∆mz values.")
            end
            mz2 = convert(Float64, method.arg[1] + method.arg[2])
            push!(xic, add_ion_current(scan.mz, scan.int, mz1, mz2))
        elseif method isa MZ
            mz1 = convert(Float64, method.arg[1])
            mz2 = convert(Float64, method.arg[2])
            push!(xic, add_ion_current(scan.mz, scan.int, mz1, mz2))
        else  # TIC
            push!(xic, scan.tic)
        end
    end

    if !isempty(xic)
        return Chromatogram(xrt, xic, maximum(xic))
    else
        ErrorException("No matching spectra.")
    end
end





"""
    average(filename::String, arguments::FilterType...; stats::Bool=true)
Returns the average mass spectrum as an [`MSscans`](@ref) container, along with the sample standard deviation of the intensities when `stats=true` (default). The data may be filtered using [`FilterType`](@ref) arguments.

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

# Examples
```julia-repl
julia> spectrum = average("test.mzXML")
MassJ.MSscans([1, 2, 3, ...

julia> spectrum = average("test.mzML", MassJ.Level(2))
MassJ.MSscans([2, 3], ...

julia> spectrum = average("test.mzXML", MassJ.Precursor(1255.5), MassJ.RT([1, 60]))
MassJ.MSscans([2, 5, 8, 11], ...
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

    elseif ext in ("mzml", "mgf", "msp", "imzml")
        scans = load(filename)
        return average(scans, arguments...; stats=stats)

    else
        ErrorException("File format not supported.")
    end
end


"""
    average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
Returns the average mass spectrum container (MSscans) along with the sample standard deviation of the intensities with stats=true (default) for all the mass spectra within the Array of mass spectrum container MSscan.. The data may be filtered by level, precursor mass, activation methods, etc, using the arguments MassJ.Level(N), MassJ.Precursor(mz), MassJ.Activation_Method("method"), or any combination of these arguments.
# Examples
```julia-repl
julia> spectrum = average("test.mzxml")
MassJ.MSscans([1, 2, 3 ....
julia> spectrum = average("test.mzxml", MassJ.Level(1) )
MassJ.MSscans([1, 4, 7, 10,
julia> spectrum = average("test.mzxml", MassJ.Precursor(1255.5) )
MassJ.MSscans([2, 5, 8, 11, ...
julia> spectrum = average("test.mzxml", MassJ.Activation_Method("PQD") )
MassJ.MSscans([3, 6, 9, 12, 15,
julia> spectrum = average("test.mzxml", MassJ.Activation_Method("PQD"), MassJ.Polarity("+"), MassJ.RT([10,20]))
MassJ.MSscans([9, 12, 15, 18], ...
```
"""
function average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
    pred = compose_predicates(scans, arguments)

    result = nothing
    count = 0
    for scan in scans
        pred(scan) || continue
        count += 1
        if count == 1
            result = scan
        elseif count == 2
            result = stats ? avg(result, scan) : result + scan
        else
            result = stats ? avg(result, scan) : result + scan
        end
    end

    if count >= 2
        return stats ? standard_deviation(result, count) : result / count
    elseif count == 1
        return result
    else
        ErrorException("No matching spectra.")
    end
end

