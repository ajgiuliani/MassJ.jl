
"""
Deconvolution submodule.
"""

# User Interface.
# ---------------

export deconv


"""
    deconv(scan::MScontainer, method::Charges; FWHM::Real=-1, R::Real=-1, shape::Symbol=:autoguess, maxiter::Int=250, tol::Real=1e-06)
Returns the result of the deconvolution of an MSscan(s) object using the UniDec algorithm. The output is an object of the same type as the input.  The deconv function takes an MSscan or MSscans object and a Charges method as mandatory inputs. Optional keyword arguments are:
    - FWHM: full width at half maximum of the peaks.
    - R: mass resolving power at m/z 500.
    - if neither R nor FWHM are provided, the peak width is obtained from the base peak FWHM of the mass spectrum.
    - shape: peak shape function: :gauss, :lorentz, :voigt or :autoguess (default).
    - maxiter: maximum number of iterations. Default value is 250.
    - tol: convergence criterion of the iterative procedure. Default is 1e-06.

# Examples
```julia-repl
julia> deconv_data = deconv(scans, Charges(adduct="H", range=(1,10)))
MassJ.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
```julia-repl
julia> deconv_data = deconv(scans, Charges(adduct="H", range=(5,15), width=2), R=5000)
MassJ.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function deconv(scan::MScontainer, method::Charges; FWHM::Real = -1, R::Real = -1, shape::Symbol =:autoguess, maxiter::Int = 250, tol::Real = 1e-06)

    charges = method.range
    zwidth = method.width
    adduct = method.adduct

    # checking for linearity of the mz axis & resampling if not
    if is_evenly_spaced(scan.mz)
        mz = scan.mz
        int = scan.int
    else
        mz, int = resampling(scan.mz, scan.int)
    end

    if length(adduct) == 0
        ma = 0.0
    else
        ma = masses(formula(adduct))["Monoisotopic"]  # monoisotopic mass of the adduct
    end

    model = _resolve_shape(scan, shape)
    FWHM = _resolve_FWHM(scan, model, R, FWHM)

    maxi = maximum(int)

    y = _deconv(mz, int .* (100.0 / maxi), charges, ma, zwidth, FWHM, model, maxiter, tol)

    maximum(charges) + 1 - minimum(charges) > 10 ? N = 10 : N = maximum(charges) + 1 - minimum(charges)

    a = minimum((mz .- ma) .* minimum(charges))
    b = maximum((mz .- ma) .* maximum(charges))
    deconv_mz = range(a, stop = b, length = length(mz) * N)

    deconv_int = zeros(maximum(charges), length(mz) * N)

    for i = minimum(charges):maximum(charges)
        temp_mz = (mz .- ma) .* i
        extrap = LinearInterpolation(temp_mz, y[i,:], extrapolation_bc = Line())
        deconv_int[i,:] = [extrap(x) for x in deconv_mz]
    end
    proj = vec(sum(deconv_int[minimum(charges):maximum(charges), :], dims=1))
    basePeakIntensity = maximum(proj)
    basePeakMz = deconv_mz[num2pnt(proj, basePeakIntensity)]
    tic = sum(proj)

    if scan isa MSscan
        return MSscan(scan.num, scan.rt, tic, collect(deconv_mz), proj .* (100.0 / basePeakIntensity), scan.level, basePeakMz, 100.0, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    elseif scan isa MSscans
        return MSscans(scan.num, scan.rt, tic, collect(deconv_mz), proj .* (100.0 / basePeakIntensity), scan.level, basePeakMz, 100.0, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, [0.0])
    end

end


"""
    _resolve_shape(scan::MScontainer, shape::Symbol)
Resolve the peak shape model from the input symbol. If :autoguess, determine the best shape from the base peak.
"""
function _resolve_shape(scan::MScontainer, shape::Symbol)
    if shape == :autoguess
        model = autoguess_shape(scan)
        println("Found $model peak shape.")
        return model
    elseif shape == :gauss
        return MassJ.gauss
    elseif shape == :lorentz
        return MassJ.lorentz
    elseif shape == :voigt
        return MassJ.voigt
    else
        error("Unknown shape. Please choose between :gauss, :lorentz, :voigt or :autoguess (default)")
    end
end


"""
    _resolve_FWHM(scan::MScontainer, model::Function, R::Real, FWHM::Real)
Resolve the peak FWHM from either the resolving power R, an explicit FWHM, or auto-detection.
"""
function _resolve_FWHM(scan::MScontainer, model::Function, R::Real, FWHM::Real)
    if R == -1 && FWHM == -1
        FWHM = autoguess_FWHM(scan, model)
        st = @sprintf("Found peak FWHM to be %.2f Da.", FWHM)
        println(st)
        return FWHM
    elseif FWHM == -1 && R != -1
        return 500.0 / R
    elseif R != -1 && FWHM != -1
        error("Please provide either the resolving power (R @ m/z 500) or the FWHM, but not both.")
    else
        return FWHM
    end
end


"""
    _deconv(mz::AbstractVector, int::AbstractVector, charges::Tuple{Int,Int}, ma::Float64, zwidth::Int, FWHM::Real, shape::Function, maxiter::Int, tol::Real)
Private function performing the UniDec deconvolution on the int and mz input arrays. It returns a matrix with the charge states as rows and the m/z points as columns.
"""
@noinline function _deconv(mz::AbstractVector, int::AbstractVector, charges::Tuple{Int,Int}, ma::Float64, zwidth::Int, FWHM::Real, shape::Function, maxiter::Int, tol::Real)

    zrange = maximum(charges)
    nmz = length(mz)

    # initialization: each charge row starts as a copy of the input spectrum
    f = zeros(zrange, nmz)
    for i = 1:zrange
        f[i,:] = int
    end
    fprev = similar(f)

    peak_shape = get_peak_shape(shape, FWHM, mz)

    p = ProgressUnknown("Iterations:")
    t0 = time()
    iter = 0
    while iter < maxiter
        copyto!(fprev, f)

        s = chargefilter(mz, int, f, charges, ma, zwidth)

        c = project_N_convolve(s, peak_shape)

        @simd for e = 1:size(f, 1)
            @inbounds f[e,:] = @. s[e,:] * (int / c)
        end

        SSD = sum(abs2, f .- fprev)
        SSD /= sum(abs2, f)
        if SSD < tol
            break
        end
        iter += 1
        ProgressMeter.next!(p)
    end
    st = @sprintf("Convergence in %.2f s after %i iterations.", time() - t0, iter)
    println(st)

    return f
end


"""
    project_N_convolve(s::AbstractMatrix, g::AbstractVector)
Project the matrix along the charge dimension and convolve the result with the peak shape function.
"""
@noinline function project_N_convolve(s::AbstractMatrix, g::AbstractVector)
    proj = vec(sum(s, dims = 1))
    direct_convolution(proj, g, origin = div(length(g), 2))
end

function new_mass(m_z::Real, e::Int, i::Int, ma::Float64)
    return ((m_z * e) + (ma * i)) / (e + i)
end


"""
    chargefilter(mz::AbstractVector, int::AbstractVector, f::AbstractMatrix, charges::Tuple{Int,Int}, ma::Float64, zw::Int)
Charge state averaging filter. For each charge state, computes the geometric mean of the signal across neighboring charge states using interpolation.
"""
@noinline function chargefilter(mz::AbstractVector, int::AbstractVector, f::AbstractMatrix, charges::Tuple{Int,Int}, ma::Float64, zw::Int)
    zrange = maximum(charges)
    s = zeros(zrange, length(mz))
    b = maximum(mz)
    a = minimum(mz)

    thres = 0.0

    Threads.@threads for e = minimum(charges):maximum(charges)
        ln = [el > thres ? log(el) : -12.0 for el in f[e,:]]
        for i = -zw:zw
            if i != 0
                if (e + i) >= 1 && (e + i) <= size(f, 1)
                    for k in eachindex(mz)
                        new_mz = new_mass(mz[k], e, i, ma)
                        if new_mz >= a && new_mz <= b
                            i0, i1 = MassJ.bracketnearest(mz, new_mz)
                            new_int = lerp(f[e+i, i0], f[e+i, i1], (new_mz - mz[i0]) / (mz[i1] - mz[i0]))
                            if new_int > thres
                                @inbounds ln[k] += log(new_int)
                            else
                                @inbounds ln[k] += -12.0
                            end
                        end
                    end
                end
            end
        end
        @inbounds s[e,:] = @. exp((1 / (2 * zw + 1)) * ln)
    end
    return s
end


"""
    get_peak_shape(shape::Function, FWHM::Real, mz::AbstractVector)
Calculate the peak shape kernel from the input shape function, FWHM and the m/z scale. Returns a normalized vector.
"""
function get_peak_shape(shape::Function, FWHM::Real, mz::AbstractVector)
    stp = mz[2] - mz[1]

    nb_pnts = Int(div(2.0 * FWHM, stp)) + 1

    g_mz = range(-1.0 * FWHM, step = stp, length = nb_pnts)

    p = [FWHM, 0.0, 1.0, 0.0]
    @inbounds g = [shape(i, p) for i in g_mz]

    g = g ./ sum(g)
end


function resampling(X::AbstractArray, Y::AbstractArray)
    newX = range(X[1], stop = X[end], length = length(X))
    return collect(newX), [LinearInterpolation(X, Y)(x) for x in newX]
end

function is_evenly_spaced(X::AbstractArray)
    step = X[2] - X[1]
    for i = 1:length(X)-1
        if abs((X[i+1] - X[i]) - step) / abs(step) >= 1e-12
            return false
        end
    end
    return true
end

function figure_of_merit(h::AbstractArray, c::AbstractArray)
    h_bar = mean(h)
    return 1 - sum((h .- c).^2) / sum((h .- h_bar).^2)
end


function autoguess_shape(scan::MScontainer)
    w = roughguess_FWHM(scan)
    m0  = scan.basePeakMz - 4.0 * w
    m1 = scan.basePeakMz + 4.0 * w
    i0 = num2pnt(scan.mz, m0)
    i1 = num2pnt(scan.mz, m1)
    p0 = [w, scan.basePeakMz, scan.basePeakIntensity, 0.0]
    y = MassJ.gauss(scan.mz[i0:i1], p0)
    var_gauss = Statistics.cor(scan.int[i0:i1], y)
    y = MassJ.lorentz(scan.mz[i0:i1], p0)
    var_lor = Statistics.cor(scan.int[i0:i1], y)
    if var_gauss >= var_lor
        model = MassJ.gauss
    else
        model = MassJ.lorentz
    end
end

function roughguess_FWHM(scan::MScontainer)
    # get highest peak
    maxi = maximum(scan.int)
    imax = num2pnt(scan.int, maxi)
    i = 1
    while scan.int[imax-i] > maxi / 2.0
        i += 1
    end
    j = 1
    while scan.int[imax+j] > maxi / 2.0
        j += 1
    end
    mz1 = lin_interp(scan.int[imax-i], scan.int[imax-i+1], scan.mz[imax-i], scan.mz[imax-i+1], maxi / 2.0)
    mz2 = lin_interp(scan.int[imax+j-1], scan.int[imax+j], scan.mz[imax+j-1], scan.mz[imax+j], maxi / 2.0)
    return (mz2 - mz1)
end

function autoguess_FWHM(scan::MScontainer, model::Function)
    w = roughguess_FWHM(scan)
    m0  = scan.basePeakMz - 4.0 * w
    m1 = scan.basePeakMz + 4.0 * w
    i0 = num2pnt(scan.mz, m0)
    i1 = num2pnt(scan.mz, m1)
    p0 = [w, scan.basePeakMz, scan.basePeakIntensity, 0.0]

    fit = curve_fit(model, scan.mz[i0:i1], scan.int[i0:i1], p0)
    return fit.param[1] * 1.1
end
