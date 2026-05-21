# Exporting data

The [`save`](@ref) function writes MassJ data back to standard
mass-spectrometry file formats. The output format is selected automatically
from the file extension, mirroring how [`load`](@ref) reads:

| Extension | Format    | Writer                |
|-----------|-----------|-----------------------|
| `.mzML`   | mzML      | [`MassJ.save_mzml`](@ref)  |
| `.mzXML`  | mzXML     | [`MassJ.save_mzxml`](@ref) |

`save` accepts any of the spectrum containers — a single
[`MassJ.MSscan`](@ref), an averaged [`MassJ.MSscans`](@ref), or a
`Vector{MSscan}` produced by [`load`](@ref).

## Basic usage

```julia-repl
julia> scans = load("input.mzML");

julia> save(scans, "output.mzML")               # round-trip
"output.mzML"

julia> save(scans, "output.mzXML")              # convert format
"output.mzXML"

julia> save(scans[1], "single_scan.mzML")       # one scan
"single_scan.mzML"
```

## Options

Both writers accept the same two keywords:

* `precision` — `64` (default, `Float64` arrays) or `32` (`Float32`, halves the
  binary payload size at the cost of ≈7 digits of precision).
* `compress` — `true` (default) to zlib-compress the binary arrays. Disable for
  faster writing of small files, at the cost of larger output.

```julia
save(scans, "out.mzML";  precision = 32)            # smaller, lossy
save(scans, "out.mzML";  compress  = false)         # plain base64, no zlib
save(scans, "out.mzXML"; precision = 32, compress = false)
```

## Round-trip fidelity

Calling `load` on the file produced by `save` recovers the same spectrum data.
The following fields are guaranteed to round-trip exactly when `precision = 64`:

- `mz`, `int` (m/z and intensity arrays)
- `level`, `polarity`, `rt`, `tic`
- `basePeakMz`, `basePeakIntensity`
- `precursor`, `chargeState`, `activationMethod`, `collisionEnergy`
- `spectrumType` (mzML only — mzXML does not encode this)

Format-level metadata (`fileDescription`, `instrumentConfiguration`,
`dataProcessing`, the full `cvList`, etc.) is emitted in a *minimal-but-valid*
form. Downstream tools that depend on rich provenance fields should expect
those slots to be blank or contain a marker indicating the file was written by
MassJ.

### Type-symmetric round-trip

`load` returns a value of the same type that was passed to `save`:

| `typeof(save argument)` | `typeof(load result)`   |
|-------------------------|-------------------------|
| `MSscan`                | `MSscan`                |
| `MSscans`               | `MSscans`               |
| `Vector{MSscan}`        | `Vector{MSscan}`        |
| `Vector{MSscans}`       | `Vector{MSscans}`       |
| (file not saved by MassJ) | `Vector{MSscan}`      |

The single-value paths (`MSscan`, `MSscans`) tag the spectrum with a MassJ
`userParam` so `load` knows to unwrap it from the surrounding container. For
files that did not pass through `save`, `load` keeps its long-standing
`Vector{MSscan}` return type.

Saving a `Vector{MSscans}` writes one spectrum per element, each carrying its
own variance array and history but without the scalar marker. `load` then
returns a `Vector{MSscans}` of the same length:

```julia
vec = [average(scans1), average(scans2), average(scans3)]   # Vector{MSscans}
save(vec, "batch.mzML")
back = load("batch.mzML")                                     # Vector{MSscans}, length 3
```

```julia-repl
julia> v       = load("input.mzML");

julia> mean_s  = average(v);                  # ::MSscans, .s populated

julia> save(mean_s, "averaged.mzML");

julia> back    = load("averaged.mzML");

julia> typeof(back) == typeof(mean_s)         # bit-symmetric
true

julia> back.s == mean_s.s                     # variance preserved exactly
true

julia> back.num == mean_s.num                 # full history preserved
true
```

### What's encoded for an averaged spectrum

Saving an `MSscans` writes more than just the mean intensity:

* The per-m/z variance (`.s`) travels as an extra `<binaryDataArray>`
  (mzML) or a second `<peaks pairOrder="variance">` child (mzXML).
* The ten vector-valued provenance fields — `num`, `rt`, `level`,
  `precursor`, `polarity`, `activationMethod`, `collisionEnergy`,
  `chargeState`, `driftTime`, `compensationVoltage` — are stored as pipe-
  separated strings in MassJ-specific `userParam`s / custom attributes.

Other tools that read these files see only the standard mean intensity and
silently ignore the MassJ extensions. The file remains a valid mzML / mzXML
for them.

The byte order differs between the two formats:

- **mzML** uses little-endian, with separate `<binaryDataArray>` elements for
  m/z and intensity.
- **mzXML** uses big-endian ("network") byte order with a single `<peaks>`
  blob containing interleaved `(m/z, intensity, m/z, intensity, …)` pairs.

Both details are handled transparently by `save` and `load`.

## Use cases

* **Format conversion** — read an mzXML produced by an older instrument
  pipeline, save as mzML for a modern tool:
  ```julia
  save(load("legacy.mzXML"), "modern.mzML")
  ```
* **Sharing processed spectra** — average, smooth, and centroid a series of
  scans, then save the result for a collaborator:
  ```julia
  avg = average("input.mzML")
  smoothed = smooth(avg)
  save(smoothed, "processed.mzML")
  ```
* **Reducing file size** — for a tabulated peak list where exact intensity
  values are not critical:
  ```julia
  save(scans, "compact.mzML"; precision = 32)
  ```
