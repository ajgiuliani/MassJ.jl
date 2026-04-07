"""
Utility functions to work on MScontainer data types.

"""
    
# User Interface.
# --------------

export avg, num2pnt



# Overloaded Base function
# ------------------------

"""
    /(a::MSscan, N::Real)
Divide the intensity and the tic data of a MSscan by a number.
```julia-repl
julia> scans[1] / 1.0e2
MSj.MSscan(1, 0.1384, 50819.5, [140. ....
```
"""
function /(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end


"""
    /(a::MSscans, N::Real)
Divide the intensity, tic and variance of a MSscans by a number.
```julia-repl
julia> a / 1.0e2
MSj.MSscans(1, 0.1384, 50819.5, [140. ....
```
"""
function /(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s / N)
end

"""
    *(a::MSscan, N::Real)
Multiply the intensity and the tic data of a MSscan by a number.
```julia-repl
julia> scans[1] * 1.0e2
MSj.MSscan(1, 0.1384, 50819.5, [140. ....
```
"""
function *(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end

"""
    *(a::MSscans, N::Real)
Multiply the intensity, tic and variance of a MSscans by a number.
```julia-repl
julia> a * 1.0e2
MSj.MSscans(1, 0.1384, 50819.5, [140. ....
```
"""
function *(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s * N)
end

"""
    *(N::Real, a::MScontainer)
Commutation of multiplication of number with MScontainer.
"""
function *(N::Real, a::MScontainer)
    return a * N
end


"""
    *(a::MScontainer, b::MScontainer)
Multiplication of mass spectra elementwise.
```julia-repl
julia> a * b
MSj.MSscans([2, 5], [0.7307, 4.344
```
"""
function *(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic * b.tic
    eps = 1e-5
    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = Interpolations.LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int .* [extrap(x) for x in a.mz ]
        mz = a.mz
      
    elseif length(a.mz) < length(b.mz)
        extrap = Interpolations.LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] .* b.int
        mz = b.mz
            
    elseif length(a.mz) == length(b.mz)
        mz = a.mz
        int = a.int .* b.int
    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]

    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end


"""
    -(a::MScontainer, b::MScontainer)
Substraction of mass spectra elementwise. Negative scan num refers the 'b' MScontainer.
```julia-repl
julia> a - b
MSj.MSscans([1, 4], [0.1384, 3.7578, -0.1384, -3.7578]...
```
"""
function -(a::MScontainer, b::MScontainer)
    #num = vcat(a.num, -b.num)
    num = vcat(a.num)
    rt  = vcat(a.rt, -b.rt)
    tic = a.tic - b.tic

    s    = Vector{Float64}(undef,0)
    int  = Vector{Float64}(undef,0)
    mz   = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = Interpolations.LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int - [extrap(x) for x in a.mz ]
        mz = a.mz
        if a isa MSscans
            if b isa MSscans
                extrap2 = Interpolations.LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = a.s + [extrap2(x) for x in a.mz]
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                extrap2 = Interpolations.LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = [extrap2(x) for x in a.mz]
            end            
        end

    elseif length(a.mz) < length(b.mz)
        extrap = Interpolations.LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] - b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = Interpolations.LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            if b isa MSscans
                s = [extrap2(x) for x in b.mz] + b.s
            else
                s = [extrap2(x) for x in b.mz]
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end            
        end

    elseif length(a.mz) == length(b.mz)
        mz = a.mz
        int = a.int - b.int
        if a isa MSscans
            if b isa MSscans
                s = a.s + b.s
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end
        end         

    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]
    
    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end


"""
    +(a::MScontainer, b::MScontainer)
Addition of mass spectra elementwise.
```julia-repl
julia> scans[1] - scans[2]
MSj.MSscans([1, 2], [0.1384, 0.7307]
```
"""
function +(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic + b.tic

    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = Interpolations.LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int + [extrap(x) for x in a.mz ]
        mz = a.mz
        if a isa MSscans
            if b isa MSscans
                extrap2 = Interpolations.LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = a.s + [extrap2(x) for x in a.mz]
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                extrap2 = Interpolations.LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = [extrap2(x) for x in a.mz]
            end            
        end
      
    elseif length(a.mz) < length(b.mz)
        extrap = Interpolations.LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] + b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = Interpolations.LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            if b isa MSscans
                s = [extrap2(x) for x in b.mz] + b.s
            else
                s = [extrap2(x) for x in b.mz]
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end            
        end
            
    elseif length(a.mz) == length(b.mz)
        mz = a.mz
        int = a.int + b.int
        if a isa MSscans
            if b isa MSscans
                s = a.s + b.s
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end
        end         
    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]

    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end



### Average & Variance calculation using an incremental Welford algorithm

"""
    avg(a::MScontainer, b::MScontainer)
Returns the average of the input mass spectra and compute the variance using an incremental Welford algorithm.
```julia-repl
julia> MSj.avg(scans[1], scans[4])
MSj.MSscans([1, 4], [0.1384, 3.7578], ....
```
"""
function avg(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic + b.tic

    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
#    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        # interpolating b.int to the a.mz values and adding the result to a.int
        extrap = Interpolations.LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
#        int = a.int + [extrap(x) for x in a.mz]
        mz = a.mz
        if a isa MSscans
            m = a.int + ([extrap(x) for x in a.mz] - a.int) / (length(a.num)+1 )
            s = a.s + ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
#            s = (a.s .* a.s) + ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
        elseif a isa MSscan
            m = a.int + ([extrap(x) for x in a.mz] - a.int) / 2
            s = ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
        end
    elseif length(a.mz) < length(b.mz)
        # interpoling a.int to the b.mz values and adding the result to b.int
        extrap = Interpolations.LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
#        int = [extrap(x) for x in b.mz] + b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = Interpolations.LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            m = [extrap(x) for x in b.mz]  + (b.int -[extrap(x) for x in b.mz] ) / (length(a.num)+1 )
            s = [extrap2(x) for x in b.mz] + (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
#            s = ([extrap2(x) for x in b.mz] .* [extrap2(x) for x in b.mz]) + (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
        elseif a isa MSscan
            m = [extrap(x) for x in b.mz] + (b.int - [extrap(x) for x in b.mz] ) / 2
            s = (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
        end
    elseif length(a.mz) == length(b.mz)
        #    else
        mz = a.mz
#        int = a.int + b.int
        if a isa MSscans
            m = a.int + (b.int - a.int) / (length(a.num)+1 )
            s = a.s  + ( (b.int - m) .* (b.int - a.int) )
#            s = (a.s .* a.s)  + ( (b.int - m) .* (b.int - a.int) )
        elseif a isa MSscan
            m = a.int + (b.int - a.int) / 2
            s = (b.int - m) .* (b.int - a.int)
        end
    end

    basePeakIntensity = maximum(m)
    basePeakMz = mz[num2pnt(m, basePeakIntensity)]


    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat(a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, m, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end

function standard_deviation(a::MSscans, N::Int)
    return MSscans(a.num, a.rt, a.tic, a.mz, a.int, a.level, a.basePeakMz, a.basePeakIntensity, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, map(sqrt, (a.s / (N -1)) ) )
end



### Utility functions

"""
    add_ion_current(x::AbstractArray, y::AbstractArray, a::Real, b::Real)
Returns sum the ion current (int) within the m/z range defined by the a and b input values.
"""
function add_ion_current(x::AbstractArray, y::AbstractArray, a::Real, b::Real)
    ia = num2pnt(x, a)
    ib = num2pnt(x, b)
    return sum( y[ia:ib] )
end


"""
    num2pnt(x::AbstractArray, val::Real)
General purpose utility function used to retrieve the index of an array for which the value is closest to the input.
"""
@noinline  function num2pnt(x::AbstractArray, val::Real)
    ibest = 1
    dxbest = abs(x[ibest] - val)
    for i in eachindex(x)
        dx = abs(x[i] - val)
        if dx < dxbest
            dxbest = dx
            ibest = i
        end
    end
    return ibest
end


"""
    findnearest(X::AbstractArray, val::Real)
Returns the index of an ordered array for which the value is closest to the input.
"""
function findnearest(X::AbstractArray, val::Real)
    if val > X[end]
        return length(X)
    elseif val < X[1]
        return 1
    end    
    lo = 1
    hi = length(X)
    while lo <= hi
        mid = div(hi + lo, 2)
        if val < X[mid]
            hi = mid - 1
        elseif val > X[mid]
            lo = mid + 1
        else
            return mid
        end
    end
    return X[lo] - val < val - X[hi] ? lo : hi
    
end

"""
    bracketnearest(X::AbstractArray, val::Real)
Returns a tuple with the two indices for which the X values bracket the input value. Requires an ordered list. 
"""
function bracketnearest(X::AbstractArray, val::Real)
    if val > X[end]
        return length(X), NaN
    elseif val < X[1]
        return NaN, 1
    end    
    lo = 1
    hi = length(X)
    while lo <= hi
        mid = div(hi + lo, 2)
        if val < X[mid]
            hi = mid - 1
        elseif val > X[mid]
            lo = mid + 1
        else
            return (mid, mid)
        end
    end
    return lo < hi ? (lo,hi) : (hi,lo)    
end


"""
     lin_interp(x0::Real, x1::Real, y0::Real, y1::Real, x::Real)
Linear interpolation between the (x0,y0) and (x1, y1). Returns the y value at x.
"""
function lin_interp(x0::Real, x1::Real, y0::Real, y1::Real, x::Real)
    return (y0 * (x1 - x) + y1 * (x -x0)) / (x1 - x0)
end


"""
     lerp(v0::Real, v1::Real, t::Real)
Linear interpolation between the (v0,v0) for t in the [0,1] range. Ref: https://en.wikipedia.org/wiki/Linear_interpolation

"""
function lerp(v0::Real, v1::Real, t::Real)
    return (1 - t) * v0 + t * v1
end


"""
    savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
Savitzky-Golay filter removes high frequency noise from data. Parameters:
    int::AbstractArray
    order::Int   order of the polynomial
    window::Int  length of the window, has to be an odd number
    deriv::Int   the order of the derivative to be computed. Default = 0 leads to smoothing only.
"""
function savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
    if window % 2 != 1
        return ErrorException("Window has to be an odd number.")
    elseif window < 1
        return ErrorException("window has to be a positive number.")
    elseif window < order + 2
        return ErrorException("window is too small for the order.")
    end
    order_range = range(1, length=(order+1))
    half_window = Int( (window-1) / 2 )
    
    b = zeros(window, order+1)

    for i = 0:order
        b[:,i+1] = [x for x = -half_window:half_window].^(i)
    end
    
    m = b * LinearAlgebra.pinv(b' * b)
    coefs = m[:,deriv + 1] * factorial(deriv)
    yfirst = int[1]*ones(half_window)
    ylast = int[end]*ones(half_window)
    pad = vcat(yfirst, int, ylast)
    
    y = convolve(coefs[end:-1:1], pad)[1 + 2 * half_window : end]
   # y = conv(coefs[end:-1:1], pad)[2 * half_window + 1 : end - 2 * half_window]
    m
    return y
end


"""
    convolve(a::AbstractArray, b::AbstractArray)
Convolve arrays a and b using the Fourier transform algorithm.
"""
function convolve(a::AbstractArray, b::AbstractArray)
    na = length(a)
    nb = length(b)
    if nb > na
        u=b
        v=a            # smallest array need to be padded
    else
        u=a
        v=b            # smallest array need to be padded
    end
    if na != nb
        pad = zeros(abs(nb - na))
        padded = vcat(v, pad)
    else
        padded = v
    end
    
    real(FFTViews.ifft( FFTViews.fft(padded) .* FFTViews.fft(u)))
end


"""
    function direct_convolution(s::AbstractArray, r::AbstractArray; mode::String="constant", cval::Real=0.0, origin::Int=0)
Convolve arrays s (signal) and r (response) using the direct algorithm. Supports constant filling mode.
"""
@noinline function direct_convolution(s::AbstractArray, r::AbstractArray; mode::String="constant", cval::Real=0.0, origin::Int=0)
    ns = length(s)        # signal
    nr = length(r)        # response function
    
    if mode == "constant"
        pad = fill(cval, nr)
    else 
        error("Filling mode unsupported")
    end
        
    padded = vcat(pad, s, pad)
    #c = zeros(length(padded))
    c = zeros(ns + nr + nr)
    
     @simd for i=nr:nr+ns
        for j=1:nr
            @inbounds c[i] += padded[i+nr-j-origin] * r[nr-j+1]
        end
    end
    return c[nr:end-(nr+1)] ./= sum(r)
end


"""
    gauss(x::Real, p::AbstractArray)
Gaussian shape function.
"""
function gauss(x::Real, p::AbstractArray)
    # Gaussian shape function
    # FWHM             = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + p[3] * exp(- ( (x-p[2])^2/(2.0 * (p[1] / 2.3548)^2) ) )

    return p[4] + p[3] * exp(- ( (x-p[2])^2/(2.0 * (p[1] / 2.3548)^2) ) )
end


"""
    gauss(x::AbstractArray, p::AbstractArray)
Gaussian shape function.
"""
function gauss(x::AbstractArray, p::AbstractArray)
    # Gaussian shape function
    # FWHM             = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + p[3] * exp(- ( (x-p[2])^2/(2.0 * (p[1] / 2.3548)^2) ) )

    return @.  p[4] + p[3] * exp(- ( (x-p[2])^2/(2.0 * (p[1] / 2.3548)^2) ) )
end

"""
    lorentz(x::Real, p::AbstractArray)
Cauchy-Lorentz shape function.
"""
function lorentz(x::Real, p::AbstractArray)
    # Lorentzian shape function
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + p[3]  / (1.0 + ((x - p[2])/(p[1]/2.0))^2 )
    
    return p[4] + p[3]  / (1.0 + ((x - p[2])/(p[1]/2.0))^2 )
    
end

"""
    lorentz(x::AbstractArray, p::AbstractArray)
Cauchy-Lorentz shape function.
"""
function lorentz(x::AbstractArray, p::AbstractArray)
    # Lorentzian shape function
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + p[3]  / (1.0 + ((x - p[2])/(p[1]/2.0))^2 )
    
    return @. p[4] + p[3]  / (1.0 + ((x - p[2])/(p[1]/2.0))^2 )
    
end

"""
    voigt(x::AbstractArray, p::AbstractArray)
Pseudo-Voigt profile function used by the TBPD method.
"""
function voigt(x::AbstractArray, p::AbstractArray)
    # pseudo-voigt profile
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    
    γg = p[1] / (2.0 * sqrt(log(2.0)))
    γl = p[1] / 2.0
    γ = (γg^5 + 2.69269 * γg^4 * γl + 2.42843 * γg^3 * γl^2 + 4.47163 * γg^2 * γl^3 + 0.07842 * γg * γl^4 + γl^5)^(1/5)

    η = 1.36603 *(γl / γ) - 0.47719 * (γl / γ)^2 + 0.11116 * (γl / γ)^3
    η = 1.36603 * 0.5 - 0.47719 * (0.5)^2 + 0.11116 * (0.5)^3
    

    L(x, Γ,x0) = (1/ (π *  Γ)) / ((x-x0)^2 + Γ^2)
    G(x,Γ,x0) = exp( -( (x-x0)^2) /  Γ^2 ) / (Γ * sqrt(π))
   #return  p[4] + 4.0 * p[3]  * (  (1 - η) * G(x,γ,p[2]) ) #
   #return  p[4] +  p[3]  * ( (1 - η) * G(x,γ,p[2]) ) 
   return @. p[4] +  p[3]  * ( η * L(x,γ,p[2]) + (1 - η) * G(x,γ,p[2]) ) 
end


"""
    voigt(x::Real, p::AbstractArray)
Pseudo-Voigt profile function used by the TBPD method.
"""
function voigt(x::Real, p::AbstractArray)
    # pseudo-voigt profile
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    
    γg = p[1] / (2.0 * sqrt(log(2.0)))
    γl = p[1] / 2.0
    γ = (γg^5 + 2.69269 * γg^4 * γl + 2.42843 * γg^3 * γl^2 + 4.47163 * γg^2 * γl^3 + 0.07842 * γg * γl^4 + γl^5)^(1/5)

    η = 1.36603 *(γl / γ) - 0.47719 * (γl / γ)^2 + 0.11116 * (γl / γ)^3
    η = 1.36603 * 0.5 - 0.47719 * (0.5)^2 + 0.11116 * (0.5)^3
    

    L(x, Γ,x0) = (1/ (π *  Γ)) / ((x-x0)^2 + Γ^2)
    G(x,Γ,x0) = exp( -( (x-x0)^2) /  Γ^2 ) / (Γ * sqrt(π))
   #return  p[4] + 4.0 * p[3]  * (  (1 - η) * G(x,γ,p[2]) ) #
   #return  p[4] +  p[3]  * ( (1 - η) * G(x,γ,p[2]) ) 
   return  p[4] +  p[3]  * ( η * L(x,γ,p[2]) + (1 - η) * G(x,γ,p[2]) ) 
end



"""
    morpholaplace(input::AbstractArray, region::Int)
Performs morphological Laplacian of the input array, as defined by the addition of the dilatation and the erosion of the input array.
"""
morpholaplace(input::AbstractArray, region::Int) = dilatation(input, region) + erosion(input, region)


"""
    morphogradient(input::AbstractArray, region::Int)
Performs morphological gradient of the input array, defined by the difference between the dilatation and the erosion of the input array.
"""
morphogradient(input::AbstractArray, region::Int) = dilatation(input, region) - erosion(input, region)


"""
    tophat(input::AbstractArray, region::Int)
Performs the Top Hat of the input Array, defined by the difference between the input and its morphological opening.
"""
tophat(input::AbstractArray, region::Int) = input - opening(input, region)


"""
    bottomhat(input::AbstractArray, region::Int)
Performs the Bottom Hat of the input Array, defined by the difference between the morphological closing of the input and the input.
"""
bottomhat(input::AbstractArray, region::Int) = closing(input, region) - input


"""
    opening(input::AbstractArray, region::Int)
Performs the morphological opening of the input Array, which is the dilatation of the erosion of the input
"""
opening(input::AbstractArray, region::Int) = dilatation(erosion(input, region), region)


"""
    closing(input::AbstractArray, region::Int)
Performs the morphological closing of the input Array, which is defined as the erosion of the dilatation of the input.
"""
closing(input::AbstractArray, region::Int) = erosion(dilatation(input, region), region)



"""
    erosion(input::AbstractArray, region::Int)
Performs the morphological erosion of the input, which is the minimum-filtering over the structuring element region.
"""
erosion(input::AbstractArray, region::Int) = extremefilt(input, minimum, region)



"""
    dilatation(input::AbstractArray, region::Int)
Performs the morphological dilatation of the input, which is the maximum-filtering over the structuring element region
"""
dilatation(input::AbstractArray, region::Int) = extremefilt(input, maximum, region)


"""
    extremefilt(input::AbstractArray, minmax::Function, region::Int)
Return the erosion or the dilation of the input over the region, which the size of the structuring element.
"""
function extremefilt(input::AbstractArray, minmax::Function, region::Int)
    output = deepcopy(input)
    if region == 1
        half_dim = 1
    else
        half_dim = Int(region / 2.0)
    end
    
    for i =1:length(input)
        if i > half_dim && i < length(input) - half_dim
            @views output[i] = minmax(input[i-half_dim : i+half_dim ])
        elseif i <= half_dim
            @views output[i] = minmax(input[1 : i+half_dim])
        elseif i >= length(input) - half_dim
            @views output[i] = minmax(input[i-half_dim : end])
        end
    end
    return output
end



