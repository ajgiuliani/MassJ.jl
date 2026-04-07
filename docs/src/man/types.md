# Data types
The main data type of the package is the abstract type [`MSj.MScontainer`](@ref).

Mass spectrometry scans are stored in the following structures, inspired from the mzXML format, which is a subtype of [`MSj.MScontainer`](@ref). 
```julia
struct MSscan <: MScontainer
    num::Int                          # num
    rt::Float64                       # retentionTime
    tic::Float64                      # totIonCurrent
    mz::Vector{Float64}               # m/z
    int::Vector{Float64}              # intensity
    level::Int                        # msLevel
    basePeakMz::Float64               # basePeakMz
    basePeakIntensity::Float64        # basePeakIntensity
    precursor::Float64                # precursorMz
    polarity::String                  # polarity
    activationMethod::String          # activationMethod
    collisionEnergy::Float64          # collisionEnergy
end
```

Another subtype, [`MSj.Chromatogram`](@ref), is used to store the retention time, the ionic current and the maximum value of the ion current.

```julia
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}               # araay of retention times
    ic::Vector{Float64}               # array of ion current
    maxic::Float64                    # maximum ion current (used in plotting normalization)
end
```

Combination of mass spectra requires another subtype of [`MSj.MScontainer`](@ref) called [`MSj.MSscans`](@ref) (notice the ending s).

```julia
struct MSscans  <: MScontainer
    num::Vector{Int}                  # num
    rt::Vector{Float64}               # retentionTime
    tic::Float64                      # totIonCurrent
    mz::Vector{Float64}               # m/z
    int::Vector{Float64}              # intensity
    level::Vector{Int}                # msLevel
    basePeakMz::Float64               # basePeakMz
    basePeakIntensity::Float64        # basePeakIntensity
    precursor::Vector{Float64}        # precursorMz
    polarity::Vector{String}          # polarity
    activationMethod::Vector{String}  # activationMethod
    collisionEnergy::Vector{Float64}  # collisionEnergy
    s::Vector{Float64}                # variance
end
```
The [`MSj.MSscans`](@ref) structure is very similar to the [`MSj.MSscan`](@ref) one, except that the fields `num`, `rt`, `precursor`, `polarity`, `activationMethod` and `collisionEnergy` are vectors. The idea is to keep track of the *history* of the operations that have led to this result. For example, if a `MSscans` element is the result of the addition of two individual scans such as scans[1] + scans[2], then the `num`field of resulting `MSscans` is an array [1, 2]. The same applies to the retention time, precursor m/z, polarity, activation method and collision energy fields.

## Deconvolution method types

The deconvolution functions use dedicated method types to dispatch to the appropriate algorithm. These types are subtypes of [`MSj.MethodType`](@ref).

[`MSj.UniDec`](@ref) is a marker type for the UniDec deconvolution algorithm.

[`MSj.Charges`](@ref) specifies charge deconvolution parameters:
```julia
@with_kw struct Charges <: MethodType
    adduct::String                # adduct ion formula (e.g. "H", "Na")
    range::Tuple{Int,Int}         # charge state range (min, max)
    width::Int = 1                # charge state filter width
end
```

[`MSj.Masses`](@ref) specifies mass deconvolution parameters:
```julia
@with_kw struct Masses <: MethodType
    adduct::String                # adduct ion formula
    range::Tuple{Int,Int}         # mass range
    width::Int = 1                # mass filter width
end
```
