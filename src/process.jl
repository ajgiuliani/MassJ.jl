"""
Processing functions submodule. 
"""

# User Interface.
# ---------------

export smooth, centroid, baseline_correction


"""
    smooth(scan::MScontainer; method::MethodType=SG(5, 9))
Smooth the intensity of the input data and returns a similar structure.
# Examples
```julia-repl
julia> smoothed_data = MSj.smooth(scans)
MSj.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
    if method isa MSj.SG
        return savitzky_golay_filtering(scan, method.order, method.window, method.derivative)
    end  
end

"""
    smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
Smooth the intensity of the input data and returns a similar structure.
# Examples
```julia-repl
julia> scans = load("filename")
julia> smoothed_data = MSj.smooth(scans)
6-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6 .....
```
"""
function smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
    if method isa MSj.SG
        sm_scans = Vector{MSscan}(undef, 0)
        for el in scans
            push!(sm_scans, savitzky_golay_filtering(el, method.order, method.window, method.derivative))
        end
        return sm_scans
    end  
end



"""
    savitzky_golay_filtering(scan::MSj.MScontainer, order::Int, window::Int, deriv::Int)
Savitzky-Golay filtering to remove the high frequency noise of int data within the MSscan(s) container.
"""
function savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
    y = savitzky_golay(scan.int, order, window, deriv)
    
    basePeakIntensity = ceil(maximum(y))
    basePeakIndex = num2pnt(y, basePeakIntensity)
    basePeakMz = scan.mz[basePeakIndex]
    
    if scan isa MSscan
        return MSscan(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    elseif scan isa MSscans
        return MSscans(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
    end
end
 

"""
    centroid(scan::MScontainer; method::MethodType=MethodType=SNRA(1., 100) )
Peak picking algorithm taking a MSscan or MSscans object as input and returning an object of the same type containing the detected peaks.  Available algorithm are : Signal to Noise Ratio (SNR) and Template Based Peak Detection (TBPD). Default method is Signal to Noise Ratio Analysis (SNRA), with default threshold = 1.0 and region = 100.
# Examples
```julia-repl
julia> centroid(scans)
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
    if method isa TBPD
        ∆mz = 500.0 / method.resolution       # according to mz / ∆mz  = R, we take the value @ m/z 500
        if method.shape == :gauss
            return tbpd(scan, gauss, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :lorentz
            return tbpd(scan, lorentz, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :voigt
            return tbpd(scan, voigt, ∆mz, convert(Float64,method.threshold))
        else
            ErrorException("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
        end

    elseif method isa SNRA
        return snra(scan, method.threshold, method.region)
        
#    elseif method isa CWT()
#        return cwt(scan, method.threshold)
#    else
#        ErrorException("Unsupported method.")
    end
    
end

"""
    centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100) )
Peak picking algorithm taking an array of MSscan as input and returning an object of the same type containing the detected peaks. Available algorithm are : Signal to Noise Ratio (SNR) and Template Based Peak Detection (TBPD). Default method is Signal to Noise Ratio Analysis (SNRA), with default threshold = 1.0 and region = 100.
# Examples
```julia-repl
julia> reduced_data = centroid(scans)
6-element Array{MSj.MSscan,1}:
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100) )
    cent_scans = Vector{MSscan}(undef,0)
    if method isa TBPD
        for el in scans
            ∆mz = 500.0 / method.resolution       # according to ∆mz / mz  = R, we take the value @ m/z 500
            if method.shape == :gauss
                push!(cent_scans, tbpd(el, gauss, ∆mz, convert(Float64,method.threshold)))
            elseif method.shape == :lorentz
                push!(cent_scans, tbpd(el, lorentz, ∆mz, convert(Float64,method.threshold)))
            elseif method.shape == :voigt
                push!(cent_scans, tbpd(el, voigt, ∆mz, convert(Float64,method.threshold)))
            else
                ErrorException("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
            end
        end
        return cent_scans
    elseif method isa SNRA
        for el in scans
            push!(cent_scans, snra(el, method.threshold, method.region))
        end
        return cent_scans
#    elseif method isa CWT()
#        return cwt(scan, method.threshold)
#    else
#        ErrorException("Unsupported method.")
    end
    
end

"""
    snra(scan::MScontainer, thres::Real, region::Int)
Signal to Noise Ratio Analysis method returning the m/z and intensity of the peaks detected.
"""
function snra(scan::MScontainer, thres::Real, region::Int)
    # declaration
    peaks_mz = Vector{Float64}(undef,0)
    peaks_int = Vector{Float64}(undef,0)
    peaks_s = Vector{Float64}(undef,0)
    
    noise  = MSj.opening(scan.int, region);
    SNR = scan.int ./ noise;
        
    maxi = maximum(scan.int)    
    for i = 2:length(scan.mz)
        if scan.int[i] >=  maxi * thres / 100. 
            if SNR[i] > SNR[i-1] && SNR[i] > SNR[i+1]
                push!(peaks_int, (scan.int[i] - noise[i]))
                push!(peaks_mz, scan.mz[i])
                if scan isa MSscans
                    push!(peaks_s, scan.s[i])
                end
            end
        end
    end
    if size(peaks_int, 1) == 0
        basePeakIntensity = NaN
        basePeakMz = NaN
    else
        basePeakIntensity = maximum(peaks_int)
        basePeakMz = peaks_mz[ num2pnt(peaks_int, basePeakIntensity) ]
    end
    
    if scan isa MSscans
        return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
    elseif scan isa MSscan
        return MSscan(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    end

end

    
"""
    tbpd(scan::MSj.MScontainer, shape::Symbol,  R::Real, thres::Real)
Template based peak detection algorithm returning the m/z and intensity of the peaks detected.
"""
function tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)   #template based peak detection
    box = num2pnt(scan.mz, scan.mz[1]+0.4) - 1        # taking a box of 0.5 width m/z
    correlation = zeros(length(scan.mz))
    maxi = maximum(scan.int)
    val = 0.0
    for i = 1:1:length(scan.mz)-box
        level = scan.int[i]
        if level >=  maxi * thres / 100. 
            bkg = 0.0
            p0 = [∆mz, scan.mz[i], level, bkg]
            ydata = [model(el, p0) for el in scan.mz[i:i+box]]
            val = Statistics.cor(scan.int[i:i+box], ydata)
        else
            val = 0.0
        end
        if val >= 0.62
            correlation[i] = val
        else
            correlation[i] = 0.0
        end
    end

    peaks_mz = Vector{Float64}(undef,0)
    peaks_int = Vector{Float64}(undef,0)
    peaks_s = Vector{Float64}(undef,0)

    diff_prev = 0.0
    diff      = 0.0
    
    for i = 2:length(correlation)
        if correlation[i] > correlation[i-1] && correlation[i] > correlation[i+1]
            max_value = maximum( scan.int[i:i+2] )
            max_index = num2pnt(scan.int, max_value)
            push!(peaks_mz, scan.mz[max_index])
            push!(peaks_int, scan.int[max_index])
            if scan isa MSscans
                push!(peaks_s, scan.s[max_index])
            end
        end
    end    
    basePeakIntensity = maximum(peaks_int)
    basePeakMz = peaks_mz[ num2pnt(peaks_int, basePeakIntensity) ]
    
    if scan isa MSscans
        return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
    elseif scan isa MSscan
        return MSscan(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    end
end


    
"""
function cwt(scan::MScontainer)                                         # Continuous Wavelet Transform 
    error("Wavelets not implemented")
end
"""


"""
    baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
Baseline correction takes a MSscan or MSscans object as input and returns an object of the same type as the input with the mass spectra corrected for their base line. Available algorithms are Top Hat, Locally Weighted Error Sum of Squares regression (LOESS) and Iterative Polynomial Smoothing Algorithm (IPSA). The default method is TopHat with a structuring element width of 100 points.
# Examples
```julia-repl
julia> reduced_data = baseline_correction(scan)
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> reduced_data = baseline_correction(scans, method = MSj.LOESS(1))
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> reduced_data = baseline_correction(scans, method = MSj.IPSA(51,100))
6-element Array{MSj.MSscan,1}:
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
    if method isa TopHat
        return tophat_filter(scan, method.region)
    elseif method isa LOESS
        return loess(scan, method.iter)
    elseif method isa IPSA
        return ipsa(scan, method.width, method.maxiter)
    end
end


"""
    baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
Baseline correction takes a Vector of MSscan as input and returns a Vector of MSscan with the mass spectra corrected for their base line. Available algorithms are Top Hat, Locally Weighted Error Sum of Squares regression (LOESS) and Iterative Polynomial Smoothing Algorithm (IPSA). The default method is TopHat with a structuring element width of 100 points.
# Examples
```julia-repl
julia> reduced_data = baseline_correction(scans)
6-element Array{MSj.MSscan,1}:
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> reduced_data = baseline_correction(scans, method = MSj.LOESS(1))
6-element Array{MSj.MSscan,1}:
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> reduced_data = baseline_correction(scans, method = MSj.IPSA(51,100))
6-element Array{MSj.MSscan,1}:
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
    bl_scans = Vector{MSscan}(undef,0)
    for el in scans
        if method isa TopHat
            push!(bl_scans, tophat_filter(el, method.region))            
        elseif method isa LOESS
            push!(bl_scans, loess(el, method.iter))
        elseif method isa IPSA
            push!(bl_scans, ipsa(el, method.width, method.maxiter))
        end       
    end
    return bl_scans
end




"""
    loess(scan::MScontainer, iter::Int )
Method  taking a MSscan or MSscans object as input and returning an object of the same type with the mass spectra without their base line, using the LOESS (Locally Weighted Error Sum of Squares regression).
"""
function loess(scan::MScontainer, iter::Int )
    n = length(scan.mz) 
    r = Int(ceil( n / 2 ))
    h = [sort(abs.(scan.mz .- scan.mz[i]))[r] for i=1:n ]
    w = clamp.(abs.( ( scan.mz .- transpose(scan.mz)) ./ h), 0.0, 1.0)
    w = (1 .- w.^3).^3
    baseline = zeros(n)
    res = zeros(n)
    delta = ones(n)
    for j=1:iter 
        for i=1:n
            weight = delta .* w[:,i]
            b = [sum(weight .* scan.int), sum(weight .* (scan.int .* scan.mz))]
            A = [sum(weight), sum(weight .* scan.mz), 
                 sum(weight .* scan.mz), sum(weight .* scan.mz.^2) ]
            A = reshape(A, 2, 2)
            beta = LinearAlgebra.pinv(A) * b
            baseline[i] = beta[1] + beta[2] * scan.mz[i]
        end
        res = scan.int - baseline
        s = Statistics.median(abs.(res))
        delta = clamp.(res ./ (6.0 .* s), -1, 1)
        delta = (1 .- delta.^2).^2
    end 
    TIC = sum(res)
    basePeakIntensity = maximum(res)
    basePeakMz = scan.mz[num2pnt(scan.int,basePeakIntensity)]
    if scan isa MSscans
        return MSscans(scan.num, scan.rt, TIC, scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
    elseif scan isa MSscan
        return MSscan(scan.num, scan.rt, TIC, scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    end
end


"""
    ipsa(scan::MScontainer, width::Real, maxiter::Int)
Method  taking a MSscan or MSscans object as input and returning an object of the same type with the mass spectra without their base line, using the iterative polynomial smoothing algorithm (IPSA) baseline correction.
"""
function ipsa(scan::MScontainer, width::Real, maxiter::Int)
    if iseven(width) 
        width -= 1
    end
    #step 1
    eps = 1e-07
    input = zeros( length(scan.int) )
    res = zeros( length(scan.int) )
    #step 2
    bkg = savitzky_golay(scan.int, 0, width,0)
    bkg_old = zeros(length(scan.int))
    res = scan.int - bkg
    # step 3
    eratio_old = 0.0
    # step 4
    counter = 1
    while true

        for i = 1:length(input)
            if scan.int[i] < bkg[i]
                input[i] = scan.int[i]
            else
                input[i] = bkg[i]
            end
        end
        bkg = savitzky_golay(input, 0, width,0)
        res = scan.int - bkg ;
        eratio = norm(bkg - bkg_old) / norm(bkg)

        if (abs(eratio - eratio_old)) < eps
            break
        elseif counter > maxiter
            break
        end
        counter += 1
        eratio_old = eratio
        bkg_old = bkg
    end
    
    basePeakIntensity = ceil(maximum(res))
    basePeakIndex = num2pnt(res, basePeakIntensity)
    basePeakMz = scan.mz[basePeakIndex]

    if scan isa MSscan
        return MSscan(scan.num, scan.rt, scan.tic, scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    elseif scan isa MSscans

        return MSscans(scan.num, scan.rt, scan.tic, scan.mz, res, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
    end

end




"""
    tophat_filter(scan::MScontainer, region::Int)
Method taking a MScontainer object as input and returning an object of the same type with the mass spectra without their base line, using the TopHat filtering algorithm.
"""
function tophat_filter(scan::MScontainer, region::Int )
    TIC = sum( tophat(scan.int, region) )
    basePeakIntensity = maximum(tophat(scan.int, region))
    basePeakMz = scan.mz[num2pnt(scan.int,basePeakIntensity)]
    if scan isa MSscans
        return MSscans(scan.num, scan.rt, TIC, scan.mz, tophat(scan.int, region), scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
    elseif scan isa MSscan
        return MSscan(scan.num, scan.rt, TIC, scan.mz, tophat(scan.int, region), scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    end
    return 
end


