"""
Submodule with types and structures used to stored the data and dispatch to the right methods.
"""

# Markers used by the export / import round-trip to identify MassJ-specific
# extensions in mzML and mzXML files (variance array on MSscans + scalar/vector
# distinction so the loaded value has the same type as the saved one).
const MASSJ_CONTAINER_PARAM      = "MassJ:container_type"
const MASSJ_VARIANCE_PARAM       = "MassJ:variance_array"
const MASSJ_SCALAR_PARAM         = "MassJ:saved_as_scalar"
const MASSJ_MZXML_CONTAINER_ATTR = "MassJContainer"      # mzXML scan attribute
const MASSJ_MZXML_VARIANCE_PAIR  = "variance"            # mzXML peaks pairOrder
const MASSJ_MZXML_SCALAR_ATTR    = "MassJSavedAsScalar"  # mzXML scan attribute


### Containers

"""
    abstract type MScontainer  end
Abstract type containing any imported data belongs to the MScontainer type.
"""
abstract type MScontainer  end


"""
    struct MSscans <: MScontainer

The single spectrum container in MassJ. A raw scan read from a file and a
composite spectrum produced by [`average`](@ref) or the arithmetic operators
share one type: all per-scan *provenance* fields are vectors, so a raw scan is
simply an `MSscans` whose provenance vectors have length 1, while a composite of
`N` scans has length-`N` provenance.

    struct MSscans <: MScontainer
        num::Vector{Int}                     # scan number(s)
        rt::Vector{Float64}                  # retention time(s)
        tic::Float64                         # total ion current
        mz::Vector{Float64}                  # m/z values
        int::Vector{Float64}                 # intensity values
        level::Vector{Int}                   # MS level(s)
        basePeakMz::Float64                  # base peak m/z
        basePeakIntensity::Float64           # base peak intensity
        precursor::Vector{Float64}           # precursor m/z value(s)
        polarity::Vector{String}             # polarity/polarities
        activationMethod::Vector{String}     # activation method(s)
        collisionEnergy::Vector{Float64}     # collision energy/energies
        s::Vector{Float64}                   # per-m/z variance (empty for a single scan)
        chargeState::Vector{Int}             # precursor charge state(s) (0 = unknown)
        spectrumType::Symbol                 # :centroid, :profile, or :unknown
        driftTime::Vector{Float64}           # ion mobility drift time(s) or 1/K0 (-1.0 = absent)
        compensationVoltage::Vector{Float64} # FAIMS/DMS compensation voltage(s) (0.0 = absent)
        mobilityType::Symbol                 # :DTIMS, :TIMS, :TWIMS, :FAIMS, or :none
        metadata::Dict{String,Any}           # additional format-specific metadata
    end

The variance accumulator `s` (Welford M2) is empty for a single scan and is
populated by [`avg`](@ref) once two or more spectra are combined. Use `length(s.num)`
to obtain the number of contributing scans.

Construct a single raw scan with the scalar 12- or 18-argument forms (provenance
values are scalars and are wrapped into length-1 vectors automatically); the
13- and 19-argument vector forms build a composite spectrum directly.
"""
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
    driftTime::Vector{Float64}        # ion mobility drift times or 1/K0 (-1.0 = not present)
    compensationVoltage::Vector{Float64} # FAIMS/DMS compensation voltages (0.0 = not present)
    mobilityType::Symbol              # :DTIMS, :TIMS, :TWIMS, :FAIMS, or :none
    metadata::Dict{String,Any}        # additional format-specific metadata

    # Composite spectrum — full 19-field constructor (vector provenance).
    function MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity,
                     precursor, polarity, activationMethod, collisionEnergy, s,
                     chargeState, spectrumType, driftTime, compensationVoltage,
                     mobilityType, metadata)
        new(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity,
            precursor, polarity, activationMethod, collisionEnergy, s,
            chargeState, spectrumType, driftTime, compensationVoltage,
            mobilityType, metadata)
    end

    # Composite spectrum — 13-field constructor (defaults for the extra metadata).
    function MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity,
                     precursor, polarity, activationMethod, collisionEnergy, s)
        new(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity,
            precursor, polarity, activationMethod, collisionEnergy, s,
            fill(0, length(num)), :unknown, fill(-1.0, length(num)),
            fill(0.0, length(num)), :none, Dict{String,Any}())
    end

    # Single raw scan — full 18-field scalar form (mirrors the former MSscan).
    function MSscans(num::Integer, rt::Real, tic::Real, mz, int, level::Integer,
                     basePeakMz::Real, basePeakIntensity::Real, precursor::Real,
                     polarity::AbstractString, activationMethod::AbstractString,
                     collisionEnergy::Real, chargeState::Integer, spectrumType::Symbol,
                     driftTime::Real, compensationVoltage::Real, mobilityType::Symbol,
                     metadata)
        new([Int(num)], [Float64(rt)], Float64(tic), mz, int, [Int(level)],
            Float64(basePeakMz), Float64(basePeakIntensity), [Float64(precursor)],
            [String(polarity)], [String(activationMethod)], [Float64(collisionEnergy)],
            Float64[], [Int(chargeState)], spectrumType, [Float64(driftTime)],
            [Float64(compensationVoltage)], mobilityType, metadata)
    end

    # Single raw scan — basic 12-field scalar form (mirrors the former MSscan).
    function MSscans(num::Integer, rt::Real, tic::Real, mz, int, level::Integer,
                     basePeakMz::Real, basePeakIntensity::Real, precursor::Real,
                     polarity::AbstractString, activationMethod::AbstractString,
                     collisionEnergy::Real)
        new([Int(num)], [Float64(rt)], Float64(tic), mz, int, [Int(level)],
            Float64(basePeakMz), Float64(basePeakIntensity), [Float64(precursor)],
            [String(polarity)], [String(activationMethod)], [Float64(collisionEnergy)],
            Float64[], [0], :unknown, [-1.0], [0.0], :none, Dict{String,Any}())
    end
end


"""
    struct IonCurrent <: MScontainer

A one-dimensional ion-current trace against a single separation axis. Replaces
the former `Chromatogram`, `Mobilogram`, and `Ionogram` types — the `axis` field
records what the abscissa `x` represents.

    struct IonCurrent <: MScontainer
        x::Vector{Float64}      # axis values (retention time, drift time, or compensation voltage)
        ic::Vector{Float64}     # ion current
        axis::Symbol            # :rt, :drift, or :cv
        mobilityType::Symbol    # :DTIMS/:TIMS/:TWIMS/:FAIMS for a mobility axis, else :none
    end

The peak ion current is computed on demand with [`maxic`](@ref) rather than
stored. Construct directly, or with the keyword form
`IonCurrent(x, ic; axis = :rt, mobilityType = :none)`.
"""
struct IonCurrent <: MScontainer
    x::Vector{Float64}
    ic::Vector{Float64}
    axis::Symbol
    mobilityType::Symbol
end

IonCurrent(x, ic; axis::Symbol = :rt, mobilityType::Symbol = :none) =
    IonCurrent(collect(Float64, x), collect(Float64, ic), axis, mobilityType)

"""
    maxic(c::IonCurrent)

Maximum ion current of the trace (`0.0` for an empty trace). Computed on demand;
replaces the former stored `Chromatogram.maxic` field.
"""
maxic(c::IonCurrent) = isempty(c.ic) ? 0.0 : maximum(c.ic)


"""
    struct MSrun <: AbstractVector{MSscans}
A full mzML/mzXML "run" — a vector of [`MSscans`](@ref) together with the
file-level metadata (instrument, software, source file, data processing) and
any pre-computed ion-current traces stored in the source file.

    struct MSrun <: AbstractVector{MSscans}
        scans::Vector{MSscans}            # spectrum list
        metadata::Dict{String,Any}        # file-level cvParams (instrument, …)
        chromatograms::Vector{IonCurrent} # pre-computed traces
    end

`MSrun` is a subtype of `AbstractVector{MSscans}`, so it transparently supports
the standard array interface: `length(run)`, `run[i]`, iteration, and
broadcasting. Code that previously worked on a `Vector{MSscans}` returned by
[`load`](@ref) keeps working unchanged. Slicing with a range (e.g. `run[1:5]`)
returns a plain `Vector{MSscans}` — the file-level metadata is dropped.

`run.metadata` is populated by [`load`](@ref) when reading mzML, and emitted
back by [`save`](@ref) so that the round-trip preserves instrument
configuration, software list, source file information, and so on.
"""
struct MSrun <: AbstractVector{MSscans}
    scans::Vector{MSscans}
    metadata::Dict{String,Any}
    chromatograms::Vector{IonCurrent}
end

MSrun(scans::Vector{MSscans}) = MSrun(scans, Dict{String,Any}(), IonCurrent[])

# AbstractVector interface — delegate to the underlying scans vector.
Base.size(run::MSrun)            = size(run.scans)
Base.getindex(run::MSrun, i::Int)= run.scans[i]
Base.getindex(run::MSrun, r)     = run.scans[r]
Base.setindex!(run::MSrun, v, i) = (run.scans[i] = v)
Base.IndexStyle(::Type{<:MSrun}) = IndexLinear()


"""
    abstract type AbstractPeak
Supertype for peak descriptors accepted by [`yields`](@ref). Two concrete subtypes
are provided: [`Peak`](@ref) for a fixed m/z window and [`TargetPeak`](@ref) for a
target m/z that is located in each spectrum at integration time.
"""
abstract type AbstractPeak end


"""
    struct Peak <: AbstractPeak
A fixed m/z integration window with a label. Used by [`yields`](@ref) to integrate
the same window across every spectrum in a series.

    struct Peak <: AbstractPeak
        mz1::Float64    # lower m/z bound
        mz2::Float64    # upper m/z bound
        label::String   # peak label (CSV column header / legend entry)
    end

Three constructor forms:

```julia
Peak(mz1, mz2, label)                # explicit window
Peak(mz, label; tol = 0.5)           # mz ± tol (absolute Δm/z)
Peak(mz, label; ppm = 5.0)           # mz ± mz·ppm·1e-6
```

The three-argument constructor enforces `mz1 <= mz2` by swapping if needed.
"""
struct Peak <: AbstractPeak
    mz1::Float64
    mz2::Float64
    label::String
    function Peak(mz1::Real, mz2::Real, label::AbstractString)
        a, b = Float64(mz1), Float64(mz2)
        a <= b ? new(a, b, String(label)) : new(b, a, String(label))
    end
end


"""
    struct TargetPeak <: AbstractPeak
A target m/z plus search half-width, located in each spectrum at integration time.
Unlike [`Peak`](@ref), the integration window of a `TargetPeak` may differ from
file to file: each spectrum is searched in `[mz - tol, mz + tol]` and a window is
derived from the located peak according to `method`.

    struct TargetPeak <: AbstractPeak
        mz::Float64        # target m/z
        label::String      # peak label
        tol::Float64       # search half-width (absolute Δm/z)
        method::Symbol     # :local_max, :edges, or :centroid
        edges::Float64     # threshold for :edges (fraction of peak max)
    end

`method` values:
- `:local_max` (default) — `argmax(int)` in the search window; window is
  ±`tol` around the found m/z.
- `:edges` — start from the local max, walk outward while
  `int > edges * peak_max`. The integration window is the located peak's
  half-max-by-`edges` footprint.
- `:centroid` — run [`MassJ.centroid`](@ref) on the averaged spectrum and pick
  the strongest centroid in the search window; window is ±`tol` around it. The
  centroiding method is passed to `yields` via the `centroid_method` keyword.

Construct with `tol` (absolute) or `ppm` (parts per million):
```julia
TargetPeak(110.5, "fragment_a"; tol = 0.5)                       # :local_max
TargetPeak(500.05, "precursor"; ppm = 5.0, method = :edges)
TargetPeak(195.09, "M+H";       tol = 0.2, method = :centroid)
```
"""
struct TargetPeak <: AbstractPeak
    mz::Float64
    label::String
    tol::Float64
    method::Symbol
    edges::Float64
end


"""
    struct YieldCurve <: MScontainer
Data structure holding peak yields measured across a series of spectra indexed by an
external parameter `x` (e.g. photon energy, wavelength, collision energy). Built by
[`yields`](@ref); plotted directly with `plot(yc)`.

    struct YieldCurve <: MScontainer
        x::Vector{Float64}                      # external parameter (one per file)
        xlabel::String                          # x-axis label (e.g. "energy (eV)")
        yields::Matrix{Float64}                 # nfiles × npeaks integrated intensities
        yields_err::Matrix{Float64}             # nfiles × npeaks 1-σ uncertainties
        tic::Vector{Float64}                    # per-file sum of peak integrals (raw)
        tic_err::Vector{Float64}                # per-file 1-σ on tic
        found_mz::Matrix{Float64}               # nfiles × npeaks located m/z (NaN for Peak)
        labels::Vector{String}                  # peak labels (length = npeaks)
        windows::Vector{Tuple{Float64,Float64}} # nominal (mz1, mz2) for each peak
        files::Vector{String}                   # source file paths (one per row)
        metadata::Dict{String,Any}              # records normalization steps applied
    end

`yields_err[i, j]` carries the 1-σ uncertainty on `yields[i, j]`, propagated by
trapezoidal integration from the per-m/z standard error of the averaged
spectrum (SEM = `sqrt(s / (N*(N-1)))` where `N = length(spec.num)` and `s` is
the Welford accumulator stored in `MSscans`). It is `NaN` when no error
information is available (single scan, `MSscan` input). `tic_err[i]` is the
combined 1-σ on `tic[i]`. [`normalize_tic`](@ref) and [`normalize_flux`](@ref)
propagate these uncertainties through their respective transformations.

`found_mz[i, j]` carries the actually-located m/z when peak `j` is a
[`TargetPeak`](@ref), and `NaN` when peak `j` is a fixed [`Peak`](@ref). For
[`TargetPeak`](@ref)s, `windows[j]` is the nominal search window
`(mz - tol, mz + tol)` — the per-file integration window may be narrower
(`:local_max`, `:centroid`) or wider (`:edges`).
"""
struct YieldCurve <: MScontainer
    x::Vector{Float64}
    xlabel::String
    yields::Matrix{Float64}
    yields_err::Matrix{Float64}
    tic::Vector{Float64}
    tic_err::Vector{Float64}
    found_mz::Matrix{Float64}
    labels::Vector{String}
    windows::Vector{Tuple{Float64,Float64}}
    files::Vector{String}
    metadata::Dict{String,Any}
end


### Methods
"""
    abstract type MethodType  end
Type containing all the methods used for filtering the data.
"""
abstract type MethodType  end


# chromatogram

"""
    struct BasePeak <: MethodType
Structure for multiple dispatching to retrieve base peak chromatogram.
"""
struct BasePeak <: MethodType
   #field = "base peak"
   BasePeak() = new()
end


"""
    struct TIC <: MethodType
Dispatching to retrieve total ion current chromatogram.
"""
struct TIC <: MethodType
   #field = "TIC"
   TIC() = new()
end


"""
    struct ∆MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
Structure for multiple dispatching to retrieve extracted ion current chromatogram around an m/z ± ∆mz value given by arg = [mz, ∆mz]
"""
struct ∆MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "∆mz range"
   ∆MZ(arg::argT) where{argT} = new{argT}(arg)
end


"""
    struct MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
Structure for multiple dispatching to retrieve extracted ion current chromatogram around for m/z in the range arg = [mz1, mz2].
"""
struct MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "mz range"
   MZ(arg::argT) where{argT} = new{argT}(arg)
end


"""
    struct SG{argT <: Int} <: MethodType   #Savitzky-Golay filtering
Structure for multiple dispatching to Savitzky-Golay filtering, providing the order, window size and derivative to be performed.  Defaults values are provided in functions calls.
"""
struct SG{argT <: Int} <: MethodType   #Savitzky-Golay filtering
    order::argT
    window::argT
    derivative::argT
    SG(order::argT, window::argT, derivative::argT) where{argT} = new{argT}(order, window, derivative)
end


"""
    struct TBPD{argT1 <: Symbol, argT2 <: Real, argT3 <: Real}  <: MethodType
Structure for multiple dispatching to Template Based Peak Detection centroiding, providing the shape of the template function, the resolution and threshold.  Defaults values are provided in functions calls.
"""
struct TBPD{argT1 <: Symbol, argT2 <: Real, argT3 <: Real}   <: MethodType
    shape::argT1
    resolution::argT2
    threshold::argT3
    TBPD(shape::argT1, resolution::argT2, threshold::argT3) where{argT1, argT2, argT3} = new{argT1, argT2, argT3}(shape, resolution, threshold)
end


"""
    struct SNRA{argT1 <: Real, argT2 <: Int}  <: MethodType
Structure for multiple dispatching to Signal to Noise Ratio Analysis centroiding, providing the threshold value and the size of the region.  Defaults values are provided in functions calls.
"""
struct SNRA{argT1 <: Real, argT2 <: Int}  <: MethodType
    threshold::argT1
    region::argT2
    SNRA(threshold::argT1, region::argT2) where{argT1, argT2} = new{argT1, argT2}(threshold, region)
end


"""
struct CWT{argT <: Real}  <: MethodType
    threshold::argT
    CWT(threshold::argT) where{argT} = new{argT}(threshold)
end
"""



"""
    TopHat{argT <: Int} <: MethodType
Structure for multiple dispatching to TopHat baseline correction. Region is used to specify the dimension over which this operation is performed.
"""
struct TopHat{argT <: Int} <: MethodType
    region::argT
    TopHat(region::argT) where{argT} = new{argT}(region)
end

"""
    LOESS{argT <: Int} <: MethodType
Structure for multiple dispatching to LOcally Weighted Error Sum of Squares regression (LOESS) baseline correction.
"""
struct LOESS{argT <: Int} <: MethodType
    iter::argT
    LOESS(iter::argT) where{argT} = new{argT}(iter)
end


"""
    struct IPSA{argT1 <: Int, argT2 <: Int} <: MethodType
Structure for multiple dispatching to iterative polynomial smoothing algorithm (IPSA) baseline correction.
"""
struct IPSA{argT1 <: Int, argT2 <: Int} <: MethodType
    width::argT1
    maxiter::argT2
    IPSA(width::argT1,maxiter::argT2) where{argT1, argT2} = new{argT1, argT2}(width, maxiter)
end


"""
    struct UniDec <: MethodType
Structure for multiple dispatching to UniDec deconvolution algorithm.
"""
struct UniDec  <: MethodType
    UniDec() = new()
end


"""
    struct Charges <: MethodType
Structure for multiple dispatching to UniDec charge deconvolution algorithm.
"""
@with_kw struct Charges <: MethodType
    adduct::String
    range::Tuple{Int,Int}
    width::Int = 1
end



"""
    struct Masses <: MethodType
Structure for multiple dispatching to UniDec mass deconvolution algorithm.
"""
@with_kw struct Masses <: MethodType
    adduct::String
    range::Tuple{Int,Int}
    width::Int = 1
end



    



### Filters

"""
    abstract type FilterType end
This type contains  the structures for filtering the data.
"""
abstract type FilterType end


"""
    RT{argT <: Union{Real, AbstractVector{<:Real},  AbstractVector{<:AbstractVector{<:Real}} }}
This type contains  the structures for filtering the data.
"""
struct RT{argT <: Union{Real, AbstractVector{<:Real},  AbstractVector{<:AbstractVector{<:Real}} }} <: FilterType
   arg::argT
   RT(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct IC{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
Used for multiple dispatching to Template Base Peak Detection centroiding, providing the shape of the template function, the resolution and threshold.  Defaults values are provided in functions calls.
"""
struct IC{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   IC(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct Level{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
Used to dispatch filters to MS level.
"""
struct Level{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
   arg::argT
   #field = "level"
   Level(arg::argT) where{argT} = new{argT}(arg)
end

"""
     Scan{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
Dispatch filter to scan num.
"""
struct Scan{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
   arg::argT
   #field = "num"
   Scan(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct Polarity{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
Dispatch filter to polarity.
"""
struct Polarity{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
   arg::argT
   #field = "polarity"
   Polarity(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct Activation_Method{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
Dispatch filter to activation methods
"""
struct Activation_Method{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
   arg::argT
   #field = "activationMethod"
   Activation_Method(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct Activation_Energy{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
Dispatch filter to activation energies.
"""
struct Activation_Energy{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   #field = "collisionEnergy"
   Activation_Energy(arg::argT) where{argT} = new{argT}(arg)
end


"""
    struct Precursor{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
Dispatch filter to precursor.
"""
struct Precursor{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   #field = "precursorMz"
   Precursor(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct DriftTime{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
Dispatch filter to ion mobility drift time or 1/K0 values.
"""
struct DriftTime{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   DriftTime(arg::argT) where{argT} = new{argT}(arg)
end

"""
    struct CompensationVoltage{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
Dispatch filter to FAIMS/DMS compensation voltage.
"""
struct CompensationVoltage{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   CompensationVoltage(arg::argT) where{argT} = new{argT}(arg)
end



# ----------------------------------------------------------------------------
# REPL-friendly `show` methods
#
# The default Julia `show(::IO, ::MIME"text/plain", x)` for a struct dumps
# every field — for `MSscans` that floods the REPL with arrays tens of
# thousands of elements long. These compact summaries print a one-liner with
# the most informative fields instead.
# ----------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", s::MSscans)
    n    = length(s.num)
    mzlo = isempty(s.mz) ? "—" : string(round(first(s.mz), digits = 3))
    mzhi = isempty(s.mz) ? "—" : string(round(last(s.mz),  digits = 3))
    if n == 1
        pol  = isempty(s.polarity) ? "?" : first(s.polarity)
        actm = isempty(s.activationMethod) || isempty(first(s.activationMethod)) ?
               "" : " " * first(s.activationMethod)
        ce   = first(s.collisionEnergy) > 0 ? "@" * string(first(s.collisionEnergy)) * "eV" : ""
        prec = first(s.precursor) > 0 ? "  precursor=" * string(first(s.precursor)) : ""
        print(io, "MSscans(num=", first(s.num),
                  ", MS", first(s.level), pol, actm, ce,
                  ", ", length(s.mz), " pts m/z=[", mzlo, ", ", mzhi, "]",
                  ", rt=", round(first(s.rt), digits = 4), " min, tic=", s.tic,
                  prec, ")")
    else
        print(io, "MSscans(composite of ", n, " scans",
                  ", ", length(s.mz), " pts m/z=[", mzlo, ", ", mzhi, "]",
                  ", tic=", s.tic, isempty(s.s) ? "" : ", with variance", ")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", c::IonCurrent)
    unit = c.axis === :rt ? " min" : ""
    xlo  = isempty(c.x) ? "—" : string(round(first(c.x), digits = 3))
    xhi  = isempty(c.x) ? "—" : string(round(last(c.x),  digits = 3))
    print(io, "IonCurrent(", c.axis, ", ", length(c.x), " pts x=[", xlo, ", ", xhi, "]", unit,
              ", max ic=", maxic(c), ")")
end

function Base.show(io::IO, ::MIME"text/plain", run::MSrun)
    keys_summary = isempty(run.metadata) ? "—" :
                   join(sort(collect(keys(run.metadata))), ", ")
    println(io, "MSrun(", length(run.scans), " scans, ",
                length(run.chromatograms), " chromatograms)")
    println(io, "  file_metadata keys: ", keys_summary)
    if !isempty(run.scans)
        s1, sn = first(run.scans), last(run.scans)
        println(io, "  first scan: ", sprint(show, MIME"text/plain"(), s1))
        if length(run.scans) > 1
            println(io, "  last  scan: ", sprint(show, MIME"text/plain"(), sn))
        end
    end
end
