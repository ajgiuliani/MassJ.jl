"""
Module for importing and exporting data. Dispatch to specific methods according to the file extension.
"""


# User Interface.
# ---------------

export info, load, chromatogram, mobilogram, ionogram, average



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
        error("File format not supported.")
    end
end



"""
    load(filename::String)
Loads the mass spectra from a file and returns an [`MSrun`](@ref) — a vector-like container holding one [`MSscans`](@ref) per scan, plus the file-level metadata and any pre-computed ion-current traces. The function dispatches to the appropriate reader based on the file extension.

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

# Examples
```julia-repl
julia> run = load("test.mzXML")
MSrun(6 scans, 0 chromatograms)
  file_metadata keys: —
  first scan: MSscans(num=1, MS1+, 22320 pts m/z=[140.083, 2000.0], rt=0.1384 min, tic=5.08195e6)
  last  scan: MSscans(num=6, MS3+ PQD@35.0eV, 23400 pts m/z=[50.083, 2000.0], rt=5.7689 min, tic=4.84455  precursor=902.33)

julia> run = load("test.mzML")
MSrun(3 scans, 0 chromatograms)
  file_metadata keys: data_processing, file_content, instruments, software, source_files
  first scan: MSscans(num=1, MS1+, 5 pts m/z=[100.0, 500.0], rt=0.5 min, tic=19000.0)
  last  scan: MSscans(num=3, MS2+ HCD@30.0eV, 3 pts m/z=[120.0, 220.0], rt=1.5 min, tic=2100.0  precursor=500.0)

julia> run = load("test.mgf")
MSrun(3 scans, 0 chromatograms)
  file_metadata keys: —
  first scan: MSscans(num=1, MS2+, 4 pts m/z=[110.0, 250.0], rt=0.5 min, tic=4800.0  precursor=400.0)
  last  scan: MSscans(num=3, MS2+, 4 pts m/z=[130.0, 310.0], rt=1.5 min, tic=2600.0  precursor=600.0)

julia> run = load("test.msp")
MSrun(3 scans, 0 chromatograms)
  file_metadata keys: —
  first scan: MSscans(num=1, MS2+@20.0eV, 5 pts m/z=[42.034, 195.088], rt=0.0 min, tic=178600.0  precursor=195.0877)
  last  scan: MSscans(num=3, MS1+, 4 pts m/z=[180.063, 183.074], rt=0.0 min, tic=8955.0)

julia> run = load("test.imzML")
MSrun(4 scans, 0 chromatograms)
  file_metadata keys: —
  first scan: MSscans(num=1, MS1+, 3 pts m/z=[100.0, 300.0], rt=0.0 min, tic=8000.0)
  last  scan: MSscans(num=4, MS1+, 3 pts m/z=[150.0, 350.0], rt=0.0 min, tic=8500.0)
```
"""
function load(filename::String; recursive::Bool = false, type = nothing)
    # A directory passed as a plain string is gathered like load([dir]).
    isdir(filename) && return load([filename]; recursive = recursive, type = type)
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
        error("File format not supported.")
    end
end


# Lowercase extension of `path` without the leading dot.
function _ext_of(path::AbstractString)
    e = lowercase(splitext(path)[2])
    return startswith(e, ".") ? e[2:end] : e
end

# True when `path` has a supported spectrum extension.
_is_supported_spectrum(path::AbstractString) = _ext_of(path) in _YIELDS_SUPPORTED_EXT

# Normalise the `type` selector to a list of allowed extensions. `nothing` means
# all supported formats; otherwise a Symbol/String or a collection of them, each
# validated against the supported set (so `:raw` is caught early).
function _allowed_exts(type)
    type === nothing && return _YIELDS_SUPPORTED_EXT
    syms = type isa Union{Symbol,AbstractString} ? (type,) : type
    exts = String[]
    for s in syms
        e = lowercase(String(s))
        e = startswith(e, ".") ? e[2:end] : e
        e in _YIELDS_SUPPORTED_EXT ||
            error("load: unsupported file type :$s (supported: $(join(_YIELDS_SUPPORTED_EXT, ", ")))")
        push!(exts, e)
    end
    return exts
end

# Collect spectrum-file paths from a list of files and/or directories, restricted
# to `type` (default: every supported format). Directories are listed in
# natural-sort order; with `recursive`, sub-directories are walked too. Explicit
# files of a *supported but unselected* type are silently skipped (the `type`
# selector's purpose); a genuinely unsupported explicit file is an error.
function _gather_spectrum_files(paths; recursive::Bool = false, type = nothing)
    allowed = _allowed_exts(type)
    files = String[]
    for p in paths
        if isdir(p)
            if recursive
                for (base, _, fs) in walkdir(p), fn in fs
                    full = joinpath(base, fn)
                    _ext_of(full) in allowed && push!(files, full)
                end
            else
                append!(files, list_spectra(p; type = type))
            end
        elseif isfile(p)
            ext = _ext_of(p)
            if ext in allowed
                push!(files, String(p))
            elseif ext in _YIELDS_SUPPORTED_EXT
                # supported format, but excluded by `type` — skip quietly
            else
                error("load: unsupported file extension: $p")
            end
        else
            error("load: path not found: $p")
        end
    end
    recursive && sort!(files, by = _natkey)
    return files
end

"""
    load(paths; recursive=false, type=nothing) -> Vector{MSscans}
Load every supported spectrum found across `paths` — a `Vector`, `Tuple`, or any
iterable collection of strings, each a spectrum file **or a directory** — returning
one combined `Vector{MSscans}`. (A single file or directory can also be passed as a
plain `String`: `load("run_A/")`.) Directories are scanned in natural-sort order for
supported files (mzML, mzXML, MGF, MSP, imzML, TXT), recursing into sub-directories
when `recursive = true`. File-level metadata (the `MSrun` wrapper) is dropped — the
result is a flat spectrum list.

`type` restricts loading to one or more formats when a folder holds several — a
single `Symbol`/`String` (`type = :mzml`) or a collection (`type = (:mzml, :mzxml)`);
files of other (even supported) formats are then ignored. `nothing` (default) loads
every supported format.

# Examples
```julia-repl
julia> scans = load("run_A/");                        # one folder (string)

julia> scans = load(["run_A/", "run_B/"]);            # several folders
julia> scans = load(("run_A/", "run_B/"));            # a tuple works too

julia> scans = load("mixed/"; type = :mzML);          # only mzML files in the folder
julia> scans = load(["batch/"]; recursive = true);    # walk sub-folders too
```
"""
function load(paths::AbstractVector{<:AbstractString}; recursive::Bool = false, type = nothing)
    files = _gather_spectrum_files(paths; recursive = recursive, type = type)
    isempty(files) && error("load: no supported spectra found in $(join(paths, ", ")).")
    out = MSscans[]
    for f in files
        append!(out, load(f))
    end
    return out
end

# Tuples / sets / other string iterables → delegate to the vector form.
load(paths::Union{Tuple{Vararg{AbstractString}},AbstractSet{<:AbstractString}};
     recursive::Bool = false, type = nothing) =
    load(collect(String, paths); recursive = recursive, type = type)


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
            error("Not an mzXML file.")
        end

        msRun = find_element(xroot, "msRun")
        scanCount = attribute(msRun, "scanCount")

        rt = retention_time(msRun)
        free(xdoc)
    elseif ext == "mzml"
        rt = retention_time_mzml(filename)
    elseif ext == "mgf" || ext == "msp" || ext == "imzml"
        scans = load(filename)
        rt = [s.rt[1] for s in scans]
    else
        error("File format not supported.")
    end
    return rt
end


"""
    retention_time(scans::AbstractVector{MSscans})
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
function retention_time(scans::AbstractVector{MSscans})
    rt  = Vector{Float64}(undef,0)
    for elem in scans
        push!(rt, elem.rt[1])
    end
    return rt
end



"""
    chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
Returns an [`IonCurrent`](@ref) holding the retention time (its `x` axis) and the ion current (`ic`) for all the mass spectra within the file; the maximum value is available on demand via [`maxic`](@ref). The `method` keyword controls how the ion current is computed: `TIC()` (default) for total ion current, `BasePeak()` for base peak intensity, `∆MZ([mz, Δ])` for extracted ion current around mz ± Δ, or `MZ([mz1, mz2])` for a m/z range.

The data may be filtered using [`FilterType`](@ref) arguments: `Level(N)`, `Precursor(mz)`, `Activation_Method("method")`, `Activation_Energy(ce)`, `Polarity("+")`, `Scan(n)`, `RT([t1,t2])`, `IC([min,max])`, `DriftTime(dt)`, `CompensationVoltage(cv)`, `ChargeState(z)`, `SpectrumType(:centroid)`, `MobilityType(:TIMS)`, `MetaData(key[, value])`, `InstrumentConfig("id")`. Filters compose (AND) in a single pass — see [Combining and filtering data](@ref).

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

# Examples
```julia-repl
julia> chrom = chromatogram("test.mzXML")
IonCurrent(rt, 6 pts x=[0.138, 5.769] min, max ic=5.08195e6)

julia> chrom = chromatogram("test.mzML", method = MassJ.BasePeak())
IonCurrent(rt, 3 pts x=[0.5, 1.5] min, max ic=8000.0)

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
            error("Not an mzXML file.")
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
            error("No matching spectra.")
        end

    elseif ext in ("mzml", "mgf", "msp", "imzml")
        # Load all scans and delegate to the AbstractVector{MSscans} method
        scans = load(filename)
        return chromatogram(scans, filters...; method=method)

    else
        error("File format not supported.")
    end
end

# Per-scan ion-current value for a chromatogram/mobilogram/ionogram `method`:
# TIC (default), base-peak intensity, or summed intensity over an m/z window
# (∆MZ = mz ± Δ, MZ = [mz1, mz2]).
function _ion_current(scan::MSscans, method::MethodType)
    if method isa BasePeak
        return scan.basePeakIntensity
    elseif method isa ∆MZ
        mz1 = convert(Float64, method.arg[1] - method.arg[2])
        mz1 < 0.0 && error("Bad mz ± ∆mz values.")
        mz2 = convert(Float64, method.arg[1] + method.arg[2])
        return add_ion_current(scan.mz, scan.int, mz1, mz2)
    elseif method isa MZ
        return add_ion_current(scan.mz, scan.int,
                               convert(Float64, method.arg[1]), convert(Float64, method.arg[2]))
    else  # TIC
        return scan.tic
    end
end

"""
    chromatogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
Returns the retention time and the total ion current by default for all the mass spectra within the Array of mass spectrum container `MSscans`. Alternatively, other options may be supplied such as method = MassJ.BasePeak, which returs the base peak intensity, method = MassJ.∆MZ([500,5]), which returns the ion current for the range mz = 500 ± 5, or method = MassJ.MZ([200,1000]) which return the ion current in the range from m/z 200 to m/z 1000.  The data may be filtered by ms level, precursor mass, activation methods, etc, using the arguments MassJ.Level(N), MassJ.Precursor(mz), MassJ.Activation_Method("method")...
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
function chromatogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
    pred = compose_predicates(scans, filters)

    xrt = Vector{Float64}(undef, 0)
    xic = Vector{Float64}(undef, 0)

    for scan in scans
        pred(scan) || continue
        push!(xrt, scan.rt[1])
        push!(xic, _ion_current(scan, method))
    end

    if !isempty(xic)
        return IonCurrent(xrt, xic; axis = :rt)
    else
        error("No matching spectra.")
    end
end


"""
    mobilogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
    mobilogram(filename::String, filters::FilterType...; method::MethodType=TIC())
Extract an ion-mobility trace (an [`IonCurrent`](@ref) on the `:drift` axis): the
drift time / `1/K₀` of each matching scan versus its ion current. The `method`
keyword controls the ion current exactly as for [`chromatogram`](@ref) (`TIC()`,
`BasePeak()`, `∆MZ([mz, Δ])`, `MZ([mz1, mz2])`), and [`FilterType`](@ref)
arguments compose as usual. Scans without a drift-time dimension
(`driftTime == -1.0`) are skipped; points follow scan order.

# Examples
```julia-repl
julia> m = mobilogram(scans, MassJ.Level(1), method = MassJ.∆MZ([885.5, 0.5]));
```
"""
function mobilogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
    pred = compose_predicates(scans, filters)
    x  = Float64[]
    ic = Float64[]
    mt = :none
    for scan in scans
        pred(scan) || continue
        scan.driftTime[1] == -1.0 && continue        # no mobility dimension
        push!(x, scan.driftTime[1])
        push!(ic, _ion_current(scan, method))
        mt = scan.mobilityType
    end
    isempty(ic) && error("No matching spectra with drift-time information.")
    return IonCurrent(x, ic; axis = :drift, mobilityType = mt)
end

mobilogram(filename::String, filters::FilterType...; method::MethodType=TIC()) =
    mobilogram(load(filename), filters...; method = method)


"""
    ionogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
    ionogram(filename::String, filters::FilterType...; method::MethodType=TIC())
Extract a differential-mobility (FAIMS/DMS) trace (an [`IonCurrent`](@ref) on the
`:cv` axis): the compensation voltage of each matching scan versus its ion
current. `method` and [`FilterType`](@ref) arguments behave as for
[`chromatogram`](@ref). Scans without a compensation-voltage dimension
(`compensationVoltage == 0.0`) are skipped; points follow scan order.

# Examples
```julia-repl
julia> g = ionogram(scans, method = MassJ.BasePeak());
```
"""
function ionogram(scans::AbstractVector{MSscans}, filters::FilterType...; method::MethodType=TIC())
    pred = compose_predicates(scans, filters)
    x  = Float64[]
    ic = Float64[]
    mt = :none
    for scan in scans
        pred(scan) || continue
        scan.compensationVoltage[1] == 0.0 && continue   # no CV dimension (sentinel)
        push!(x, scan.compensationVoltage[1])
        push!(ic, _ion_current(scan, method))
        mt = scan.mobilityType
    end
    isempty(ic) && error("No matching spectra with compensation-voltage information.")
    return IonCurrent(x, ic; axis = :cv, mobilityType = mt)
end

ionogram(filename::String, filters::FilterType...; method::MethodType=TIC()) =
    ionogram(load(filename), filters...; method = method)





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
function average(filename::String, arguments::FilterType...;
                 stats::Bool=true, recursive::Bool=false, type = nothing)
    # A directory passed as a plain string is grand-averaged like average([dir]).
    isdir(filename) && return average([filename], arguments...;
                                      stats = stats, recursive = recursive, type = type)
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
        if length(indices) >= 2
            return composite_spectra(filename, indices, stats)
        elseif length(indices) == 1
            return load_mzxml(filename, indices[1])
        else
            error("No matching spectra.")
        end

    elseif ext in ("mzml", "mgf", "msp", "imzml")
        scans = load(filename)
        return average(scans, arguments...; stats=stats)

    else
        error("File format not supported.")
    end
end


"""
    average(scans::AbstractVector{MSscans}, arguments::FilterType...; stats::Bool=true)
Returns the average mass spectrum container (`MSscans`) along with the sample standard deviation of the intensities with stats=true (default) for all the mass spectra within the Array of mass spectrum container `MSscans`. The data may be filtered by level, precursor mass, activation methods, etc, using the arguments MassJ.Level(N), MassJ.Precursor(mz), MassJ.Activation_Method("method"), or any combination of these arguments.
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
function average(scans::AbstractVector{MSscans}, arguments::FilterType...; stats::Bool=true)
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
        error("No matching spectra.")
    end
end


"""
    average(paths::AbstractVector{<:AbstractString}, arguments::FilterType...; stats=true, recursive=false)
Grand-average every supported spectrum found across `paths` (folders and/or
files), after optional [`FilterType`](@ref) filtering. Paths are gathered as in
[`load(::AbstractVector)`](@ref). Returns a single [`MSscans`](@ref).

# Examples
```julia-repl
julia> spec = average(["run_A/", "run_B/"], MassJ.Level(1));
```
"""
function average(paths::AbstractVector{<:AbstractString}, arguments::FilterType...;
                 stats::Bool = true, recursive::Bool = false, type = nothing)
    scans = load(paths; recursive = recursive, type = type)
    return average(scans, arguments...; stats = stats)
end

# Tuples / sets / other string iterables → delegate to the vector form.
average(paths::Union{Tuple{Vararg{AbstractString}},AbstractSet{<:AbstractString}},
        arguments::FilterType...; stats::Bool = true, recursive::Bool = false, type = nothing) =
    average(collect(String, paths), arguments...;
            stats = stats, recursive = recursive, type = type)

