# Combining and filtering data
## Average

The [`average`](@ref) returns the average of the mass spectra directly from a `Vector{MSscans}` after [Importing data](@ref) data or directly from the `filename`.

```julia-repl
julia> average("filename")
MassJ.MSscans([1, 2, 3 ....

julia> scans = load("filename")
51-element Array{MassJ.MSscans,1}:
 MassJ.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
...

julia> average(scans)
MassJ.MSscans([1, 2, 3 ....

```

Operating on files takes more time than working on `Vector{MSscans}` but may be useful to reduce the memory load.

Without any argument the [`average`](@ref) function averages the entire content of the data and the [`chromatogram`](@ref) function operates on also on the entire data.


## Filtering

The [`average`](@ref) and [`chromatogram`](@ref) functions may takes arguments to select specific fields of interest within the data and operate on them. The argument belongs to the [`MassJ.FilterType`](@ref). Their properties are listed below:

| FilterType                 | Description           | Arguments                                | Specificity            |
|----------------------------|-----------------------|------------------------------------------|------------------------|
| MassJ.Scan                   | Scan num              | Int, Vector{Int}                         | average, chromatogram |
| MassJ.Level                  | MS level              | Int, Vector{Int}                         | average, chromatogram |
| MassJ.Polarity               | Polarity              | String, Vector{String}                   | average, chromatogram |
| MassJ.Activation_Method      | Activation method     | String, Vector{String}                   | average, chromatogram |
| MassJ.Activation_Energy      | Activation energy     | Real, Vector{Real}                       | average, chromatogram |
| MassJ.Precursor              | Precursor _m/z_       | Real, Vector{Real}                       | average, chromatogram |
| MassJ.RT                     | Retention time        | Real, Vector{Real}, Vector{Vector{Real}} | average               |
| MassJ.IC                     | Ion current           | Vector{Real}                             | average               |
| MassJ.DriftTime              | Ion mobility drift time | Real, Vector{Real}                     | average, chromatogram |
| MassJ.CompensationVoltage    | FAIMS/DMS CV          | Real, Vector{Real}                       | average, chromatogram |
| MassJ.ChargeState            | Precursor charge state | Int, Vector{Int}                        | average, chromatogram, extract |
| MassJ.SpectrumType           | :centroid / :profile  | Symbol, Vector{Symbol}                   | average, chromatogram, extract |
| MassJ.MobilityType           | :DTIMS/:TIMS/:FAIMS/:none | Symbol, Vector{Symbol}               | average, chromatogram, extract |
| MassJ.MetaData               | per-spectrum cvParam  | key + value / [lo,hi] / substring / key-only | average, chromatogram, extract |
| MassJ.InstrumentConfig       | instrument configuration | String (id / cvParam name / component) | average, chromatogram, extract |

The last block is new in v2 and targets the richer per-spectrum information that the mzML
reader exposes. Every `FilterType` is converted to a predicate and **composed with `all`**, so
any combination of filters — old or new — can be passed together and is applied in a single pass.

### Filtering on charge state, spectrum type and mobility

These three filters select on dedicated [`MassJ.MSscans`](@ref) fields. Like every filter they
accept a single value or a vector of values:

```julia
extract(run, ChargeState(2))                     # only 2+ precursors
extract(run, ChargeState([2, 3]))                # 2+ or 3+
extract(run, SpectrumType(:centroid))            # only centroided spectra
extract(run, MobilityType(:TIMS))                # only TIMS scans
extract(run, MobilityType([:TIMS, :FAIMS]))      # TIMS or FAIMS
```

### Filtering on acquisition cvParams — `MetaData`

`MassJ.MetaData` filters on any key the reader stored in a spectrum's `metadata` dictionary
(`filter_string`, `ion_injection_time`, `mass_resolving_power`, `scan_window_lower`/`upper`,
isolation-window parameters, `spectrum_title`, …). The match depends on the value you supply:

```julia
MetaData("mass_resolving_power", 60000.0)        # numeric exact match
MetaData("mass_resolving_power", [50000, Inf])   # numeric range  lo ≤ x ≤ hi
MetaData("filter_string", "FTMS")                # substring match on the stored string
MetaData("ion_injection_time")                   # presence test (the key exists)
```

### Filtering by instrument configuration — `InstrumentConfig`

`MassJ.InstrumentConfig` selects spectra by the instrument configuration they were acquired
with. The query is resolved through the run-level instrument table carried by the
[`MassJ.MSrun`](@ref): it matches the configuration `id`, any of its cvParam names/accessions, or
a component type (source / analyzer / detector).

```julia
extract(run, InstrumentConfig("IC1"))         # by configuration id
extract(run, InstrumentConfig("orbitrap"))    # by a cvParam name of the configuration
extract(run, InstrumentConfig("analyzer"))    # by component type
```

(When applied to a plain `Vector{MSscans}` with no run-level table, the query is matched
directly against each spectrum's stored `instrument_configuration_ref`.)

### Combining old and new filters

New and existing filters compose freely — pass as many as you like and they are AND-combined in
a single pass. For example, to average the centroided 2+ MS² spectra acquired on the orbitrap at
high resolving power:

```julia
run = load("data.mzML")
average(run, Level(2),
             SpectrumType(:centroid),
             ChargeState(2),
             InstrumentConfig("orbitrap"),
             MetaData("mass_resolving_power", [60000, Inf]))
```



When the argument is restricted to a single value, such as `MassJ.Scan(1)`, filtering is performed on that specific value. If the argument is a vector then filtering involves all the values within the range.  Filtering on `MassJ.scan([1,10])` means that the result will be obtained for scans ranging from 1 to 10.  The same applies for all `FilterType` with the exception of `MassJ.∆MZ`, for which the first value of the vector represents the *mz* and the second value represents the spread ∆mz, so that filtering is operated for all *mz* value in the range [m/z - ∆mz , m/z + ∆mz].  The `MassJ.RT` type may take a vector or vectors as argument, such `MassJ.RT([ [1,10], [20, 30] ]).  In that case, mass spectra will be averaged in [1,10] and [20,30] range.


These filters may be combined together if necessary. For example, the input below returns the average mass spectrum for:
- the MS2 scans (level = 2), 
- precursor m/z 1255.5, 
- upon CID activation conditions
- with an activation energy of 18 
- and for retention times in the 1 to 60 s range.

```julia
average("filename", MassJ.Precursor(1255.5),
                    MassJ.Activation_Energy(18),
                    MassJ.Activation_Method("CID"),
                    MassJ.Level(2),
                    MassJ.RT( [1, 60] ),
                    )
```

Several filter types may also be combined for `chromatograms`:
```julia
chromatogram("filename", MassJ.Precursor(1255.5),
                         MassJ.Activation_Energy(18),
                         MassJ.Activation_Method("CID"),
                         MassJ.Level(2),
                         )
```

If the condition does not match any existing data, then an `ErrorException` is returned with the `"No matching spectra."` message.


The `chromatogram` function has some methods using `MassJ.MethodType` arguments:

| MethodType   | Description         | Arguments    | Remark  |
|--------------|---------------------|--------------|---------|
| MassJ.TIC      | Total ion current   | None         | Default |
| MassJ.BasePeak | Base peak intensity | None         |         |
| MassJ.MZ       | *m/z* range         | Vector{Real} |         |
| MassJ.∆MZ      | *m/z* ± ∆mz         | Vector{Real} |         |


These types control the way chromatograms are calculated: either using the total ionic current, the base peak intensity or using a *m/z* range.  The `method` argument of the `MassJ.chromatogram`function is set to MassJ.TIC() by default. This setting may be overruled by setting the method to desired value:

```julia
chromatogram("filename", method = MassJ.BasePeak())
chromatogram("filename", method = MassJ.MZ( [257, 259] ) ) 
chromatogram("filename", method = MassJ.∆MZ( [258, 1] ) ) 
```

## Extracting subsets

The [`extract`](@ref) returns a Vector of `MSscans`from either a file of from a Vector{MSscans} following a [`load`](@ref) command, which corresponds to the filter conditions. See the [Filtering](@ref) part above.

```julia
sub_set = extract("filename")                       # extracting without any conditions returns a vector identical to the output 
sub_set = extract("filename", MassJ.Level(2) )        # extract MS/MS spectra
scans = load("test.mzxml")                          # load mass spectra
sub_set = extract(scans)                            # extract a sub_set without conditions returns the original data
```
