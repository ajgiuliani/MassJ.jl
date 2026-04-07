# Data types
The main data type of the package is the abstract type [`MSj.MScontainer`](@ref).

Mass spectrometry scans are stored in the following structure, which is a subtype of [`MSj.MScontainer`](@ref).
```julia
struct MSscan <: MScontainer
    num::Int                          # scan number
    rt::Float64                       # retention time (minutes)
    tic::Float64                      # total ion current
    mz::Vector{Float64}              # m/z values
    int::Vector{Float64}             # intensity values
    level::Int                        # MS level
    basePeakMz::Float64              # base peak m/z
    basePeakIntensity::Float64       # base peak intensity
    precursor::Float64               # precursor m/z
    polarity::String                 # polarity ("+" or "-")
    activationMethod::String         # activation method (e.g. "CID", "HCD")
    collisionEnergy::Float64         # collision energy
    chargeState::Int                 # precursor charge state (0 = unknown)
    spectrumType::Symbol             # :centroid, :profile, or :unknown
    driftTime::Float64               # ion mobility drift time or 1/K0 (-1.0 = not present)
    compensationVoltage::Float64     # FAIMS/DMS compensation voltage (0.0 = not present)
    mobilityType::Symbol             # :DTIMS, :TIMS, :TWIMS, :FAIMS, or :none
    metadata::Dict{String,Any}       # additional format-specific metadata
end
```

A backward-compatible constructor accepting the original 12 fields is provided. The 6 new fields (`chargeState`, `spectrumType`, `driftTime`, `compensationVoltage`, `mobilityType`, `metadata`) default to neutral values (0, `:unknown`, -1.0, 0.0, `:none`, empty dict).

Another subtype, [`MSj.Chromatogram`](@ref), is used to store the retention time, the ionic current and the maximum value of the ion current.

```julia
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}               # array of retention times
    ic::Vector{Float64}               # array of ion current
    maxic::Float64                    # maximum ion current (used in plotting normalization)
end
```

## Ion mobility container types

Two additional container types are provided for ion mobility data:

[`MSj.Mobilogram`](@ref) stores drift time vs intensity data (analogous to a chromatogram for ion mobility):
```julia
struct Mobilogram <: MScontainer
    dt::Vector{Float64}               # drift time or 1/K0 values
    ic::Vector{Float64}               # ion current
    maxic::Float64                    # maximum ion current
    mobilityType::Symbol              # :DTIMS, :TIMS, :TWIMS, or :none
end
```

[`MSj.Ionogram`](@ref) stores compensation voltage vs intensity data (for FAIMS/DMS differential mobility):
```julia
struct Ionogram <: MScontainer
    cv::Vector{Float64}               # compensation voltage values
    ic::Vector{Float64}               # ion current
    maxic::Float64                    # maximum ion current
end
```

## Averaged spectra

Combination of mass spectra requires another subtype of [`MSj.MScontainer`](@ref) called [`MSj.MSscans`](@ref) (notice the ending s).

```julia
struct MSscans  <: MScontainer
    num::Vector{Int}                  # scan numbers
    rt::Vector{Float64}               # retention times
    tic::Float64                      # total ion current
    mz::Vector{Float64}               # m/z values
    int::Vector{Float64}              # intensity values
    level::Vector{Int}                # MS levels
    basePeakMz::Float64               # base peak m/z
    basePeakIntensity::Float64        # base peak intensity
    precursor::Vector{Float64}        # precursor m/z values
    polarity::Vector{String}          # polarities
    activationMethod::Vector{String}  # activation methods
    collisionEnergy::Vector{Float64}  # collision energies
    s::Vector{Float64}                # variance
    chargeState::Vector{Int}          # precursor charge states (0 = unknown)
    spectrumType::Symbol              # :centroid, :profile, or :unknown
    driftTime::Vector{Float64}        # ion mobility drift times (-1.0 = not present)
    compensationVoltage::Vector{Float64} # FAIMS/DMS compensation voltages (0.0 = not present)
    mobilityType::Symbol              # :DTIMS, :TIMS, :TWIMS, :FAIMS, or :none
    metadata::Dict{String,Any}        # additional format-specific metadata
end
```
The [`MSj.MSscans`](@ref) structure is similar to [`MSj.MSscan`](@ref), except that the fields `num`, `rt`, `precursor`, `polarity`, `activationMethod`, `collisionEnergy`, `chargeState`, `driftTime`, and `compensationVoltage` are vectors. This design keeps track of the *history* of the operations. For example, if an `MSscans` element is the result of the addition of two individual scans such as `scans[1] + scans[2]`, then the `num` field of the resulting `MSscans` is `[1, 2]`.

A backward-compatible constructor accepting the original 13 fields is provided. The 6 new fields default to neutral values.

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
