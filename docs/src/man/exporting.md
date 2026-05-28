# Exporting data

The [`save`](@ref) function writes MassJ data back to standard
mass-spectrometry file formats. The output format is selected automatically
from the file extension, mirroring how [`load`](@ref) reads:

| Extension | Format    | Writer                     | Notes                              |
|-----------|-----------|----------------------------|------------------------------------|
| `.mzML`   | mzML      | [`MassJ.save_mzml`](@ref)  | PSI standard; indexed by default   |
| `.mzXML`  | mzXML     | [`MassJ.save_mzxml`](@ref) | legacy, big-endian arrays          |
| `.mgf`    | MGF       | [`MassJ.save_mgf`](@ref)   | Mascot Generic Format peak lists   |
| `.msp`    | MSP       | [`MassJ.save_msp`](@ref)   | NIST spectral-library text         |
| `.txt`    | TXT       | [`MassJ.save_txt`](@ref)   | two columns; one spectrum per file |

`save` accepts any of the spectrum containers — a single or averaged
[`MassJ.MSscans`](@ref), a `Vector{MSscans}`, or the [`MassJ.MSrun`](@ref)
produced by [`load`](@ref).

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

The binary-format writers (mzML, mzXML) accept these keywords; the text writers
(mgf, msp, txt) ignore them:

* `precision` — `64` (default, `Float64` arrays) or `32` (`Float32`, halves the
  binary payload size at the cost of ≈7 digits of precision).
* `compress` — `true` (default) to zlib-compress the binary arrays. Disable for
  faster writing of small files, at the cost of larger output.
* `progress` — `true` (default) shows a `ProgressMeter` bar while writing. Set
  `false` in scripts / CI to silence it.
* `indexed` (mzML only) — `true` (default) wraps the output in `<indexedmzML>`
  with an `<indexList>` of byte offsets to each spectrum/chromatogram and an
  SHA-1 `<fileChecksum>`. Most modern proteomics tools (MaxQuant in particular)
  require this wrapper. Set `false` to emit plain `<mzML>`.

```julia
save(scans, "out.mzML";  precision = 32)            # smaller, lossy
save(scans, "out.mzML";  compress  = false)         # plain base64, no zlib
save(scans, "out.mzML";  progress  = false)         # silent
save(scans, "out.mzML";  indexed   = false)         # plain mzML (no <indexedmzML>)
save(scans, "out.mzXML"; precision = 32, compress = false)
```

## Indexed mzML output

By default the mzML writer wraps its output in the `<indexedmzML>` element
defined by the PSI mzML 1.1 indexed schema. The wrapper adds three things at
the end of the file:

1. An `<indexList>` mapping each `<spectrum>`/`<chromatogram>` id to its byte
   offset in the file, so tools that need random access can jump directly to
   a given scan without scanning the whole file.
2. An `<indexListOffset>` giving the byte offset of the index itself.
3. An SHA-1 `<fileChecksum>` computed over every byte from the start of the
   file up to and including the `<fileChecksum>` open tag.

MassJ's reader recognises both indexed and non-indexed mzML transparently
— files saved either way round-trip cleanly. Use `indexed = false` only
when the receiving tool refuses the wrapper (rare; the wrapped form is the
default for ProteoWizard's msConvert and is required by MaxQuant, among
others).

## Performance

The writers stream directly to disk one spectrum at a time, so peak RAM is
bounded by the largest single spectrum (typically tens of KB) plus the size of
the in-memory `Vector{MSscans}` you're writing — *not* by the total file size.
A 1.2 GB proteomics file that previously needed ~30 GB of RAM (and ~45 min)
now writes in a few minutes with memory bounded by the input vector you
already loaded.

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

### Round-trip return type

For the binary formats, `load` returns:

| `typeof(save argument)`       | `typeof(load result)` |
|-------------------------------|-----------------------|
| `MSscans` (a single spectrum) | `MSscans` (bare)      |
| `Vector{MSscans}` / `MSrun`   | `MSrun`               |
| (file not saved by MassJ)     | `MSrun`               |

A single spectrum saved as a scalar is tagged with a MassJ `userParam` so `load`
unwraps it back to a bare [`MassJ.MSscans`](@ref); every multi-spectrum file
loads as an [`MassJ.MSrun`](@ref) (which is an `AbstractVector{MSscans}`).

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

## Text formats (MGF / MSP / TXT)

The text writers emit the standard fields of each format and round-trip with the
corresponding reader on what that format can represent:

- **MGF** ([`MassJ.save_mgf`](@ref)) — one `BEGIN IONS … END IONS` block per
  spectrum with `TITLE`, `MSLEVEL`, `PEPMASS` (precursor), `CHARGE` (charge +
  polarity sign), `RTINSECONDS`, and `m/z intensity` peak rows.
- **MSP** ([`MassJ.save_msp`](@ref)) — one `Name: … Num Peaks: N` entry per
  spectrum with `PrecursorMZ`, `MSLEVEL`, `Ion_mode` (polarity),
  `Collision_energy`, `RT`, the peak rows, and a blank-line separator.
- **TXT** ([`MassJ.save_txt`](@ref)) — two whitespace-separated columns,
  `m/z  intensity`. A TXT file holds a **single** spectrum, so passing more than
  one spectrum is an error; pass a single [`MassJ.MSscans`](@ref) (e.g. `run[i]`
  or `average(run)`).

These formats carry only peak lists plus a handful of header fields, so
information outside that scope (variance array, full provenance vectors,
ion-mobility values, file-level metadata) is **not** written — use mzML to
preserve everything.

```julia
save(scans, "peaklist.mgf")          # MS/MS peak lists for a search engine
save(scans, "library.msp")           # build a spectral library
save(average(scans), "mean.txt")     # a single averaged spectrum, two columns
```

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
