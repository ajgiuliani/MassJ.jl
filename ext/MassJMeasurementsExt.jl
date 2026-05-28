module MassJMeasurementsExt

using MassJ
using MassJ: MSscans, YieldCurve, sem, stdev
using Measurements: measurement

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

end # module
