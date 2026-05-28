# Data types
The main data type of the package is the abstract type [`MassJ.MScontainer`](@ref).

!!! note "Changed in v2.0"
    Earlier versions had two spectrum types, `MSscan` (a raw scan) and `MSscans`
    (a composite of several scans), and three separate trace types,
    `Chromatogram`, `Mobilogram`, and `Ionogram`. Version 2.0 unifies each group
    into a single type: [`MassJ.MSscans`](@ref) for spectra and
    [`MassJ.IonCurrent`](@ref) for traces. See the migration notes in
    `CHANGELOG.md`.

## Spectra: `MSscans`

Every mass spectrum — whether a single scan read from a file or a composite
produced by [`average`](@ref) or the arithmetic operators — is stored in one
type, [`MassJ.MSscans`](@ref):

```julia
struct MSscans <: MScontainer
    num::Vector{Int}                     # scan number(s)
    rt::Vector{Float64}                  # retention time(s), minutes
    tic::Float64                         # total ion current
    mz::Vector{Float64}                  # m/z values
    int::Vector{Float64}                 # intensity values
    level::Vector{Int}                   # MS level(s)
    basePeakMz::Float64                  # base peak m/z
    basePeakIntensity::Float64           # base peak intensity
    precursor::Vector{Float64}           # precursor m/z value(s)
    polarity::Vector{String}             # polarity/polarities ("+" or "-")
    activationMethod::Vector{String}     # activation method(s) (e.g. "CID", "HCD")
    collisionEnergy::Vector{Float64}     # collision energy/energies
    s::Vector{Float64}                   # per-m/z variance (empty for a single scan)
    chargeState::Vector{Int}             # precursor charge state(s) (0 = unknown)
    spectrumType::Symbol                 # :centroid, :profile, or :unknown
    driftTime::Vector{Float64}           # ion mobility drift time(s) or 1/K0 (-1.0 = absent)
    compensationVoltage::Vector{Float64} # FAIMS/DMS compensation voltage(s) (0.0 = absent)
    mobilityType::Symbol                 # :DTIMS, :TIMS, :TWIMS, :FAIMS, or :none
    metadata::Dict{String,Any}           # additional format-specific metadata
end
```

The per-scan *provenance* fields (`num`, `rt`, `level`, `precursor`, `polarity`,
`activationMethod`, `collisionEnergy`, `chargeState`, `driftTime`,
`compensationVoltage`) are vectors. A raw scan has length-1 provenance, so
`scans[1].num` is `[1]` rather than `1`; a composite of several scans carries
one entry per contributing scan. This keeps the *history* of the operations:
the result of `scans[1] + scans[2]` has `num == [1, 2]`. Use `length(s.num)` to
get the number of contributing scans.

The variance accumulator `s` is empty for a single scan and is filled in by
[`avg`](@ref) once two or more spectra are combined (a Welford M2 accumulator).

Construct a single raw scan with the scalar 12- or 18-argument form (provenance
values are scalars and are wrapped into length-1 vectors automatically); the 13-
and 19-argument vector forms build a composite spectrum directly.

## Traces: `IonCurrent`

A one-dimensional ion-current trace against a single separation axis is stored
in [`MassJ.IonCurrent`](@ref):

```julia
struct IonCurrent <: MScontainer
    x::Vector{Float64}      # axis values (retention time, drift time, or compensation voltage)
    ic::Vector{Float64}     # ion current
    axis::Symbol            # :rt, :drift, or :cv
    mobilityType::Symbol    # :DTIMS/:TIMS/:TWIMS/:FAIMS for a mobility axis, else :none
end
```

`axis` records what the abscissa `x` represents: `:rt` for a chromatogram,
`:drift` for a mobilogram, `:cv` for a FAIMS/DMS ionogram. The peak ion current
is obtained on demand with [`MassJ.maxic`](@ref) rather than stored. The keyword
constructor `IonCurrent(x, ic; axis = :rt, mobilityType = :none)` is the usual
entry point; [`chromatogram`](@ref) returns an `:rt` trace.

## mzML run wrapper: `MSrun`

[`MassJ.MSrun`](@ref) wraps the result of [`load`](@ref) on an mzML file:
```julia
struct MSrun <: AbstractVector{MSscans}
    scans::Vector{MSscans}            # spectrum list
    metadata::Dict{String,Any}        # file-level cvParams (instrument, software, …)
    chromatograms::Vector{IonCurrent} # pre-computed traces
end
```

`MSrun` is a subtype of `AbstractVector{MSscans}`, so the standard array
interface (`length`, indexing, iteration, broadcasting) is delegated to the
underlying `scans` vector. Slicing returns a plain `Vector{MSscans}` — the
metadata is dropped because the slice no longer corresponds to a complete
run.

`metadata` is populated by [`load`](@ref) from the top-level mzML sections
(`fileDescription`, `softwareList`, `instrumentConfigurationList`,
`dataProcessingList`, `referenceableParamGroupList`) and re-emitted by
[`save`](@ref) so that the round-trip preserves instrument configuration,
source file information, processing history, and so on.

For non-mzML formats (mzXML, MGF, MSP, imzML, TXT), `load` returns
`Vector{MSscans}` directly.

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
when no error information is available (a single scan).

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
