"""
Submodule with types and structures used to stored the data and dispatch to the right methods.
"""

### Containers

"""
    abstract type MScontainer  end
Abstract type containing any imported data belongs to the MScontainer type.
"""
abstract type MScontainer  end


"""
    struct MSscan <: MScontainer
Data structure used to store individual mass spectrometry scans organized following the structure of mzXML files.

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


"""
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

"""
    struct MSscans  <: MScontainer
Data structure designed to store mass spectra obtained after filtering operation along with the history of these operation.

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

"""
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


"""
    struct Chromatogram  <: MScontainer
Data structure used to retrieve chromatography data.

    struct Chromatogram  <: MScontainer   
        rt::Vector{Float64}               # retention time
        ic::Vector{Float64}               # ion current
        maxic::Float64
    end

"""
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}
    ic::Vector{Float64}
    maxic::Float64
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

