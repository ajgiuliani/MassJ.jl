# Data types
The main data type of the package is the abstract type [`MassJ.MScontainer`](@ref).

Mass spectrometry scans are stored in the following structure, which is a subtype of [`MassJ.MScontainer`](@ref).
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

Another subtype, [`MassJ.Chromatogram`](@ref), is used to store the retention time, the ionic current and the maximum value of the ion current.

```julia
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}               # array of retention times
    ic::Vector{Float64}               # array of ion current
    maxic::Float64                    # maximum ion current (used in plotting normalization)
end
```

## Ion mobility container types

Two additional container types are provided for ion mobility data:

[`MassJ.Mobilogram`](@ref) stores drift time vs intensity data (analogous to a chromatogram for ion mobility):
```julia
struct Mobilogram <: MScontainer
    dt::Vector{Float64}               # drift time or 1/K0 values
    ic::Vector{Float64}               # ion current
    maxic::Float64                    # maximum ion current
    mobilityType::Symbol              # :DTIMS, :TIMS, :TWIMS, or :none
end
```

[`MassJ.Ionogram`](@ref) stores compensation voltage vs intensity data (for FAIMS/DMS differential mobility):
```julia
struct Ionogram <: MScontainer
    cv::Vector{Float64}               # compensation voltage values
    ic::Vector{Float64}               # ion current
    maxic::Float64                    # maximum ion current
end
```

## Averaged spectra

Combination of mass spectra requires another subtype of [`MassJ.MScontainer`](@ref) called [`MassJ.MSscans`](@ref) (notice the ending s).

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
The [`MassJ.MSscans`](@ref) structure is similar to [`MassJ.MSscan`](@ref), except that the fields `num`, `rt`, `precursor`, `polarity`, `activationMethod`, `collisionEnergy`, `chargeState`, `driftTime`, and `compensationVoltage` are vectors. This design keeps track of the *history* of the operations. For example, if an `MSscans` element is the result of the addition of two individual scans such as `scans[1] + scans[2]`, then the `num` field of the resulting `MSscans` is `[1, 2]`.

A backward-compatible constructor accepting the original 13 fields is provided. The 6 new fields default to neutral values.

## Peak and yield-curve types

[`MassJ.AbstractPeak`](@ref) is the supertype for peak descriptors accepted by
[`yields`](@ref). Two concrete subtypes are provided:

[`MassJ.Peak`](@ref) carries a fixed m/z window used identically across every
spectrum in a series:
```julia
struct Peak <: AbstractPeak
    mz1::Float64    # lower m/z bound
    mz2::Float64    # upper m/z bound
    label::String
end
```

[`MassJ.TargetPeak`](@ref) carries a target m/z and a search half-width; the
window is determined per file using one of three algorithms (`:local_max`,
`:edges`, `:centroid`):
```julia
struct TargetPeak <: AbstractPeak
    mz::Float64        # target m/z
    label::String
    tol::Float64       # search half-width (absolute Δm/z)
    method::Symbol     # :local_max, :edges, or :centroid
    edges::Float64     # threshold (fraction of max) for :edges
end
```

[`MassJ.YieldCurve`](@ref) holds the result of [`yields`](@ref):
```julia
struct YieldCurve <: MScontainer
    x::Vector{Float64}                      # external parameter, one per file
    xlabel::String                          # x-axis label (e.g. "energy (eV)")
    yields::Matrix{Float64}                 # nfiles × npeaks integrated intensities
    yields_err::Matrix{Float64}             # nfiles × npeaks 1-σ uncertainties (NaN = unknown)
    tic::Vector{Float64}                    # per-file sum of peak integrals (raw)
    tic_err::Vector{Float64}                # per-file 1-σ on tic
    found_mz::Matrix{Float64}               # nfiles × npeaks located m/z (NaN for Peak)
    labels::Vector{String}                  # peak labels
    windows::Vector{Tuple{Float64,Float64}} # nominal (mz1, mz2) for each peak
    files::Vector{String}                   # source file paths
    metadata::Dict{String,Any}              # records normalization steps applied
end
```
`yields_err` and `tic_err` carry propagated 1-σ uncertainties — see
[Uncertainties](@ref) in the energy-resolved yields manual. They are `NaN`
when no error information is available (single scan / `MSscan` input).

## Deconvolution method types

The deconvolution functions use dedicated method types to dispatch to the appropriate algorithm. These types are subtypes of [`MassJ.MethodType`](@ref).

[`MassJ.UniDec`](@ref) is a marker type for the UniDec deconvolution algorithm.

[`MassJ.Charges`](@ref) specifies charge deconvolution parameters:
```julia
@with_kw struct Charges <: MethodType
    adduct::String                # adduct ion formula (e.g. "H", "Na")
    range::Tuple{Int,Int}         # charge state range (min, max)
    width::Int = 1                # charge state filter width
end
```

[`MassJ.Masses`](@ref) specifies mass deconvolution parameters:
```julia
@with_kw struct Masses <: MethodType
    adduct::String                # adduct ion formula
    range::Tuple{Int,Int}         # mass range
    width::Int = 1                # mass filter width
end
```
