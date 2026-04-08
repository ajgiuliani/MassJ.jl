# Importing data

The [`load`](@ref) function reads a mass spectrometry file and returns a `Vector{MSscan}`. The file format is automatically determined from the extension.

Supported file formats: mzXML, mzML, MGF, MSP, imzML, TXT.

```julia-repl
julia> scans = load("test.mzXML")
6-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.1384, 5.08195e6, ...)
...

julia> scans = load("test.mzML")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.5, 19000.0, ...)
...

julia> scans = load("test.mgf")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.5, 4800.0, ...)
...

julia> scans = load("library.msp")
3-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.0, 178600.0, ...)
...

julia> scans = load("sample.imzML")
10000-element Array{MassJ.MSscan,1}:
 MassJ.MSscan(1, 0.0, 8000.0, ...)
...
```

Individual scans may be retrieved from the array the usual way:
```julia-repl
julia> scans[1]
MassJ.MSscan(1, 0.1384, 5.08195e6, ...)

julia> scans[1].mz
22320-element Array{Float64,1}:
 140.083
 ...
```

## Format-specific notes

### mzXML
The mzXML reader supports both 32-bit and 64-bit precision data, zlib compression, and nested scan hierarchies (MS1 → MS2 → MS3).

### mzML
The mzML reader supports the PSI (Proteomics Standards Initiative) standard format. It handles both `indexedmzML` and raw `mzML` files. Binary data arrays are decoded from base64 with optional zlib compression in little-endian byte order. Retention times are automatically converted to minutes regardless of the unit in the file (minutes or seconds). Ion mobility metadata (drift time, 1/K0, compensation voltage) is extracted when present.

### MGF
The MGF (Mascot Generic Format) reader loads centroided peak lists. Each `BEGIN IONS`...`END IONS` block becomes one `MSscan`. The `PEPMASS` field sets the precursor m/z, `CHARGE` sets the charge state, and `RTINSECONDS` is converted to minutes. The `TITLE` and original `SCANS` values are stored in the `metadata` dictionary.

### MSP
The MSP reader loads spectra from NIST Mass Spectral Library files used by NIST, MoNA, MassBank, and GNPS. Each entry starts with `Name:` and contains metadata fields (e.g. `Precursor_mz`, `Ion_mode`, `Collision_energy`, `Formula`, `InChIKey`) followed by `Num Peaks:` and the peak list. Both one-pair-per-line and semicolon-separated formats are supported. Rich metadata (name, formula, MW, InChIKey, CAS#, DB#, precursor type, comments) is stored in the `metadata` dictionary of each scan.

### imzML
The imzML reader loads imaging mass spectrometry data. The format consists of an XML metadata file (`.imzML`) following the mzML schema with IMS-specific CV terms for spatial coordinates, and a companion binary data file (`.ibd`) in the same directory. Both continuous and processed storage modes are supported. Spatial coordinates (x, y, and optionally z) are stored in the `metadata` dictionary as `"position_x"`, `"position_y"`, and `"position_z"`. The `info` function also reports image dimensions. Binary data supports 32-bit and 64-bit precision with optional zlib compression.

### TXT
The TXT reader loads a single spectrum from a two-column whitespace-separated file (m/z and intensity).

## Chromatograms

Chromatograms may be retrieved from a file and imported in [`MassJ.Chromatogram`](@ref):
```julia-repl
julia> chromatogram("test.mzML")
MassJ.Chromatogram([0.5, 1.0, 1.5], [19000.0, 4800.0, 2100.0], 19000.0)
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
