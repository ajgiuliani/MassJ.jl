# Information
The [`info`](@ref) public function reads the content of a file, but without loading the mass spectrometry data, and returns a `Vector{String}` containing the number of scans, scan levels, and for MS/MS data, the precursor m/z, the activation method and energy. Additional information may be gained by setting `verbose = true`.

Supported file formats: mzXML, mzML, MGF, MSP, imzML.

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

julia> info("test.mgf")
4-element Array{String,1}:
 "3 scans"
 "MS2+ 400.0"
 "MS2+ 500.0"
 "MS2+ 600.0"

julia> info("library.msp")
4-element Array{String,1}:
 "3 scans"
 "MS2+ 195.0877"
 "MS2- 179.0344"
 "MS1+"

julia> info("sample.imzML")
3-element Array{String,1}:
 "10000 spectra"
 "MS1+"
 "Image dimensions: 100 x 100 pixels"
```
