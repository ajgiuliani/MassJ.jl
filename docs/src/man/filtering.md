# Combining and filtering data
## Average

The [`average`](@ref) returns the average of the mass spectra directly from a `Vector{MSscan}` after [Importing data](@ref) data or directly from the `filename`.

```julia-repl
julia> average("filename")
MSj.MSscans([1, 2, 3 ....

julia> scans = load("filename")
51-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
...

julia> average(scans)
MSj.MSscans([1, 2, 3 ....

```

Operating on files takes more time than working on `Vector{MSscan}` but may be useful to reduce the memory load.

Without any argument the [`average`](@ref) function averages the entire content of the data and the [`chromatogram`](@ref) function operates on also on the entire data.


## Filtering

The [`average`](@ref) and [`chromatogram`](@ref) functions may takes arguments to select specific fields of interest within the data and operate on them. The argument belongs to the [`MSj.FilterType`](@ref). Their properties are listed below:

| FilterType                 | Description           | Arguments                                | Specificity            |
|----------------------------|-----------------------|------------------------------------------|------------------------|
| MSj.Scan                   | Scan num              | Int, Vector{Int}                         | average, chromatogram |
| MSj.Level                  | MS level              | Int, Vector{Int}                         | average, chromatogram |
| MSj.Polarity               | Polarity              | String, Vector{String}                   | average, chromatogram |
| MSj.Activation_Method      | Activation method     | String, Vector{String}                   | average, chromatogram |
| MSj.Activation_Energy      | Activation energy     | Real, Vector{Real}                       | average, chromatogram |
| MSj.Precursor              | Precursor _m/z_       | Real, Vector{Real}                       | average, chromatogram |
| MSj.RT                     | Retention time        | Real, Vector{Real}, Vector{Vector{Real}} | average               |
| MSj.IC                     | Ion current           | Vector{Real}                             | average               |
| MSj.DriftTime              | Ion mobility drift time | Real, Vector{Real}                     | average, chromatogram |
| MSj.CompensationVoltage    | FAIMS/DMS CV          | Real, Vector{Real}                       | average, chromatogram |



!!! note

    The filtering function goes first through all the arguments and setup an array of scan num that matches the conditions. Then it uses this array to calculate the average mass spectrum.  So this procedure needs two passes through the data, which is not very efficient. This is a point to make better in the future.
 

When the argument is restricted to a single value, such as `MSj.Scan(1)`, filtering is performed on that specific value. If the argument is a vector then filtering involves all the values within the range.  Filtering on `MSj.scan([1,10])` means that the result will be obtained for scans ranging from 1 to 10.  The same applies for all `FilterType` with the exception of `MSj.∆MZ`, for which the first value of the vector represents the *mz* and the second value represents the spread ∆mz, so that filtering is operated for all *mz* value in the range [m/z - ∆mz , m/z + ∆mz].  The `MSj.RT` type may take a vector or vectors as argument, such `MSj.RT([ [1,10], [20, 30] ]).  In that case, mass spectra will be averaged in [1,10] and [20,30] range.


These filters may be combined together if necessary. For example, the input below returns the average mass spectrum for:
- the MS2 scans (level = 2), 
- precursor m/z 1255.5, 
- upon CID activation conditions
- with an activation energy of 18 
- and for retention times in the 1 to 60 s range.

```julia
average("filename", MSj.Precursor(1255.5),
                    MSj.Activation_Energy(18),
                    MSj.Activation_Method("CID"),
                    MSj.Level(2),
                    MSj.RT( [1, 60] ),
                    )
```

Several filter types may also be combined for `chromatograms`:
```julia
chromatogram("filename", MSj.Precursor(1255.5),
                         MSj.Activation_Energy(18),
                         MSj.Activation_Method("CID"),
                         MSj.Level(2),
                         )
```

If the condition does not match any existing data, then an `ErrorException` is returned with the `"No matching spectra."` message.


The `chromatogram` function has some methods using `MSj.MethodType` arguments:

| MethodType   | Description         | Arguments    | Remark  |
|--------------|---------------------|--------------|---------|
| MSj.TIC      | Total ion current   | None         | Default |
| MSj.BasePeak | Base peak intensity | None         |         |
| MSj.MZ       | *m/z* range         | Vector{Real} |         |
| MSj.∆MZ      | *m/z* ± ∆mz         | Vector{Real} |         |


These types control the way chromatograms are calculated: either using the total ionic current, the base peak intensity or using a *m/z* range.  The `method` argument of the `MSj.chromatogram`function is set to MSj.TIC() by default. This setting may be overruled by setting the method to desired value:

```julia
chromatogram("filename", method = MSj.BasePeak())
chromatogram("filename", method = MSj.MZ( [257, 259] ) ) 
chromatogram("filename", method = MSj.∆MZ( [258, 1] ) ) 
```

## Extracting subsets

The [`extract`](@ref) returns a Vector of `MSscan`from either a file of from a Vector{MSscan} following a [`load`](@ref) command, which corresponds to the filter conditions. See the [Filtering](@ref) part above.

```julia
sub_set = extract("filename")                       # extracting without any conditions returns a vector identical to the output 
sub_set = extract("filename", MSj.Level(2) )        # extract MS/MS spectra
scans = load("test.mzxml")                          # load mass spectra
sub_set = extract(scans)                            # extract a sub_set without conditions returns the original data
```
