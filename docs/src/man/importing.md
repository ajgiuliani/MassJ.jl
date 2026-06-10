# Importing data

The [`load`](@ref) function reads a mass spectrometry file and returns an [`MassJ.MSrun`](@ref) — a vector of [`MassJ.MSscans`](@ref) bundled with any file-level metadata. The file format is automatically determined from the extension. (A TXT file, or a single spectrum saved as a scalar, comes back as a bare [`MassJ.MSscans`](@ref) instead — see [Data types](types.md).)

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

```julia-repl
julia> scans = load("test.mzXML")
6-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.1384, 5.08195e6, ...)
...

julia> scans = load("test.mzML")
3-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.5, 19000.0, ...)
...

julia> scans = load("test.mgf")
3-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.5, 4800.0, ...)
...

julia> scans = load("library.msp")
3-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.0, 178600.0, ...)
...

julia> scans = load("sample.imzML")
10000-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.0, 8000.0, ...)
...
```

Individual scans may be retrieved from the array the usual way:
```julia-repl
julia> scans[1]
MassJ.MSscans(1, 0.1384, 5.08195e6, ...)

julia> scans[1].mz
22320-element Array{Float64,1}:
 140.083
 ...
```

## Loading several files and folders

[`load`](@ref) (and [`average`](@ref)) also accept a **directory**, or a
collection of files and/or directories, and pool every supported spectrum into one
flat `Vector{MSscans}`. A single path is given as a plain string; several are
given as a `Vector`, `Tuple`, or any iterable of strings:

```julia
scans = load("run_A/")                         # one folder
scans = load(["run_A/", "run_B/"])             # several folders / files
scans = load(("run_A/", "run_B/"))             # a tuple works too
scans = load(["batch/"]; recursive = true)     # walk sub-folders as well
```

Directories are scanned in *natural-sort* order. When a folder holds more than one
supported format, restrict the load with `type` — a single `Symbol`/`String` or a
collection — and the others are ignored:

```julia
scans = load("mixed/"; type = :mzML)           # only mzML files
scans = load("mixed/"; type = (:mzml, :mzxml)) # mzML and mzXML
```

[`average`](@ref) takes the same path forms (plus its [`MassJ.FilterType`](@ref)
arguments), grand-averaging everything that is gathered:

```julia
avg = average("run_A/", MassJ.Level(1); type = :mzML)
```

## Format-specific notes

### mzXML
The mzXML reader supports both 32-bit and 64-bit precision data, zlib compression, and nested scan hierarchies (MS1 → MS2 → MS3).

### mzML
The mzML reader supports the PSI (Proteomics Standards Initiative) standard format. It handles both `indexedmzML` and raw `mzML` files. Binary data arrays are decoded from base64 with optional zlib compression in little-endian byte order. Retention times are automatically converted to minutes regardless of the unit in the file (minutes or seconds). Ion mobility metadata (drift time, 1/K0, compensation voltage) is extracted when present.

[`load`](@ref) returns a [`MassJ.MSrun`](@ref) — a wrapper that holds the
spectrum vector alongside the *file-level* metadata (instrument configuration,
software list, source file, data processing, referenceable parameter groups)
and any pre-computed chromatograms baked into the source file. The mzML reader
populates all of these; the other multi-spectrum formats leave the metadata and
chromatogram list empty. `MSrun` is a subtype of `AbstractVector{MSscans}`, so
code that treated `load`'s output as a vector keeps working unchanged:

```julia-repl
julia> run = load("input.mzML");

julia> typeof(run)
MassJ.MSrun

julia> length(run), eltype(run)
(26906, MassJ.MSscans)

julia> run[1]                                  # indexing
MSscans(…)

julia> run.metadata["instruments"][1]["id"]    # file-level metadata
"IC1"

julia> run.chromatograms                        # pre-computed TIC etc.
IonCurrent[]
```

Slicing with a range (`run[1:5]`) returns a plain `Vector{MSscans}` — the
file-level metadata is dropped because the slice no longer represents a
complete run. A TXT file (a single spectrum) and a scalar-saved single
spectrum are the exceptions: they load as a bare [`MassJ.MSscans`](@ref).

Several optional cvParams that don't have a dedicated [`MassJ.MSscans`](@ref)
field are extracted into the scan's `metadata::Dict{String,Any}`. Keys are
added only when the corresponding cvParam is present in the source file
(absent fields stay absent — no zero or `missing` placeholders):

| metadata key                          | mzML cvParam                       | type     |
|---------------------------------------|------------------------------------|----------|
| `"spectrum_title"`                    | MS:1000796 spectrum title          | `String` |
| `"lowest_observed_mz"`                | MS:1000528 lowest observed m/z     | `Float64`|
| `"highest_observed_mz"`               | MS:1000527 highest observed m/z    | `Float64`|
| `"mass_resolving_power"`              | MS:1000800 mass resolving power    | `Float64`|
| `"ion_injection_time"`                | MS:1000927 ion injection time (ms) | `Float64`|
| `"scan_window_lower"`                 | MS:1000501 scan window lower limit | `Float64`|
| `"scan_window_upper"`                 | MS:1000500 scan window upper limit | `Float64`|
| `"filter_string"`                     | MS:1000512 filter string (Thermo)  | `String` |
| `"isolation_window_target_mz"`        | MS:1000827 isolation window target | `Float64`|
| `"isolation_window_lower_offset"`     | MS:1000828 isolation window lower offset | `Float64`|
| `"isolation_window_upper_offset"`     | MS:1000829 isolation window upper offset | `Float64`|
| `"selected_ion_peak_intensity"`       | MS:1000042 peak intensity          | `Float64`|

These are commonly needed in proteomics workflows (e.g. ion injection time
for TMT/iTRAQ ratio correction, scan window bounds for QC filtering). They
also round-trip cleanly through [`save`](@ref) — saving a scan with these
keys populated emits them back as the appropriate cvParams.

```julia-repl
julia> scans = load("input.mzML");

julia> scans[1].metadata["ion_injection_time"]
9.80973815918

julia> scans[1].metadata["scan_window_upper"]
1650.0
```

### MGF
The MGF (Mascot Generic Format) reader loads centroided peak lists. Each `BEGIN IONS`...`END IONS` block becomes one `MSscans`. The `PEPMASS` field sets the precursor m/z, `CHARGE` sets the charge state, and `RTINSECONDS` is converted to minutes. The `TITLE` and original `SCANS` values are stored in the `metadata` dictionary.

### MSP
The MSP reader loads spectra from NIST Mass Spectral Library files used by NIST, MoNA, MassBank, and GNPS. Each entry starts with `Name:` and contains metadata fields (e.g. `Precursor_mz`, `Ion_mode`, `Collision_energy`, `Formula`, `InChIKey`) followed by `Num Peaks:` and the peak list. Both one-pair-per-line and semicolon-separated formats are supported. Rich metadata (name, formula, MW, InChIKey, CAS#, DB#, precursor type, comments) is stored in the `metadata` dictionary of each scan.

### imzML
The imzML reader loads imaging mass spectrometry data. The format consists of an XML metadata file (`.imzML`) following the mzML schema with IMS-specific CV terms for spatial coordinates, and a companion binary data file (`.ibd`) in the same directory. Both continuous and processed storage modes are supported. Spatial coordinates (x, y, and optionally z) are stored in the `metadata` dictionary as `"position_x"`, `"position_y"`, and `"position_z"`. The `info` function also reports image dimensions. Binary data supports 32-bit and 64-bit precision with optional zlib compression.

### TXT
The TXT reader loads a single spectrum from a two-column whitespace-separated file (m/z and intensity).

## Chromatograms

Chromatograms may be retrieved from a file and imported in [`MassJ.IonCurrent`](@ref):
```julia-repl
julia> chromatogram("test.mzML")
MassJ.IonCurrent([0.5, 1.0, 1.5], [19000.0, 4800.0, 2100.0], 19000.0)
```

## Retention time

The function [`MassJ.retention_time`](@ref) reads the retention times from a file and returns a `Vector{Float64}` in minutes.
```julia-repl
julia> MassJ.retention_time("test.mzML")
3-element Array{Float64,1}:
  0.5
  1.0
  1.5
```
