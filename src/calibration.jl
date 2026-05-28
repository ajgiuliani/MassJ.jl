"""
Mass calibration / recalibration.

Fit a polynomial that maps observed m/z onto known reference m/z (lock masses),
then apply it to a spectrum or a whole run. Polynomial calibration is linear in
its coefficients, so the fit is an ordinary least-squares solve of a Vandermonde
system rather than a nonlinear optimisation.
"""

export Calibration, calibrate

"""
    struct Calibration

A fitted m/z calibration. The correction is a polynomial in the observed m/z,

    corrected = coeffs[1] + coeffs[2]·mz + coeffs[3]·mz² + …

    struct Calibration
        model::Symbol               # :linear, :quadratic, or :poly
        coeffs::Vector{Float64}     # polynomial coefficients, low order first
        ref_mz::Vector{Float64}     # known reference (lock-mass) m/z
        obs_mz::Vector{Float64}     # measured m/z of the reference ions
        residuals_ppm::Vector{Float64} # (corrected − reference)/reference · 1e6
    end

A `Calibration` is callable: `cal(mz)` returns the corrected value for a scalar
or vector of m/z.
"""
struct Calibration
    model::Symbol
    coeffs::Vector{Float64}
    ref_mz::Vector{Float64}
    obs_mz::Vector{Float64}
    residuals_ppm::Vector{Float64}
end

# Horner evaluation of a polynomial (coefficients low order first).
function _polyval(c::AbstractVector{<:Real}, x::Real)
    v = zero(float(x))
    @inbounds for k in length(c):-1:1
        v = v * x + c[k]
    end
    return v
end

(cal::Calibration)(mz::Real) = _polyval(cal.coeffs, mz)
(cal::Calibration)(mz::AbstractVector{<:Real}) = [_polyval(cal.coeffs, x) for x in mz]

"""
    calibrate(obs_mz, ref_mz; model = :linear, degree = 0) -> Calibration

Fit a calibration mapping measured `obs_mz` onto known reference `ref_mz`
(lock masses). `model` selects the polynomial degree: `:linear` (1),
`:quadratic` (2), or `:poly` with an explicit `degree`. At least `degree + 1`
reference points are required.

# Examples
```julia-repl
julia> cal = calibrate([100.002, 200.004, 300.006], [100.0, 200.0, 300.0]);

julia> cal(250.0)            # corrected m/z
249.995...
```
"""
function calibrate(obs_mz::AbstractVector{<:Real}, ref_mz::AbstractVector{<:Real};
                   model::Symbol = :linear, degree::Integer = 0)
    length(obs_mz) == length(ref_mz) ||
        error("obs_mz and ref_mz must have equal length (got $(length(obs_mz)) and $(length(ref_mz))).")
    deg = model === :linear    ? 1 :
          model === :quadratic ? 2 :
          model === :poly      ? (degree > 0 ? Int(degree) : 2) :
          error("Unknown model $model. Use :linear, :quadratic, or :poly (with degree=).")
    n = length(obs_mz)
    n >= deg + 1 ||
        error("a degree-$deg fit needs at least $(deg + 1) reference points, got $n.")

    x = collect(Float64, obs_mz)
    y = collect(Float64, ref_mz)
    V = [xi^p for xi in x, p in 0:deg]      # Vandermonde, n × (deg+1)
    coeffs = V \ y
    fitted = V * coeffs
    residuals_ppm = (fitted .- y) ./ y .* 1e6
    return Calibration(model, coeffs, y, x, residuals_ppm)
end

"""
    calibrate(scan::MSscans, cal::Calibration) -> MSscans
    calibrate(scans::AbstractVector{MSscans}, cal::Calibration)
    calibrate(run::MSrun, cal::Calibration) -> MSrun

Apply a fitted [`Calibration`](@ref) to a spectrum, a vector of spectra, or a
whole run, correcting the m/z axis (and the base-peak m/z). Intensities,
variance, and separation/IM axes are unaffected.
"""
function calibrate(scan::MSscans, cal::Calibration)
    new_mz = cal(scan.mz)
    new_bpmz = isnan(scan.basePeakMz) ? scan.basePeakMz : cal(scan.basePeakMz)
    return MSscans(scan.num, scan.rt, scan.tic, new_mz, scan.int, scan.level,
                   new_bpmz, scan.basePeakIntensity, scan.precursor, scan.polarity,
                   scan.activationMethod, scan.collisionEnergy, scan.s, scan.chargeState,
                   scan.spectrumType, scan.driftTime, scan.compensationVoltage,
                   scan.mobilityType, scan.metadata)
end

calibrate(scans::AbstractVector{MSscans}, cal::Calibration) =
    MSscans[calibrate(s, cal) for s in scans]

calibrate(run::MSrun, cal::Calibration) =
    MSrun(MSscans[calibrate(s, cal) for s in run.scans], run.metadata, run.chromatograms)

function Base.show(io::IO, cal::Calibration)
    rms = isempty(cal.residuals_ppm) ? NaN : sqrt(mean(abs2, cal.residuals_ppm))
    print(io, "Calibration(", cal.model, ", ", length(cal.obs_mz),
          " points, RMS residual ", round(rms, digits = 3), " ppm)")
end
