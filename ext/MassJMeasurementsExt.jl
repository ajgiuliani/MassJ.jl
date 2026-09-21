module MassJMeasurementsExt

using MassJ
using MassJ: MSscans, YieldCurve, sem, stdev
using Measurements: measurement, Measurement, value, uncertainty

# intensity ± σ  (σ = standard error of the mean, or sample standard deviation)
function MassJ.measurements(spec::MSscans; kind::Symbol = :sem)
    e = kind === :sem ? sem(spec) :
        kind === :std ? stdev(spec) :
        error("measurements: `kind` must be :sem or :std (got :$kind).")
    isempty(e) &&
        error("measurements(spec): no replicate statistics — run `average(...; stats = true)` first.")
    return measurement.(spec.int, e)
end

# yield ± yields_err  (NaN where no scan-to-scan estimate is available)
MassJ.measurements(yc::YieldCurve) = measurement.(yc.yields, yc.yields_err)

# The way back: columns computed on `measurements(yc)` (any expression) are split
# into value / 1-σ and handed to the core `with_yields`.
MassJ.with_yields(yc::YieldCurve, labels, M::AbstractVecOrMat{<:Measurement}; kwargs...) =
    MassJ.with_yields(yc, labels, value.(M), uncertainty.(M); kwargs...)

end # module
