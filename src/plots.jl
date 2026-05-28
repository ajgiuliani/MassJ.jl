"""
Plotting module for MScontainer data type (MSscans and IonCurrent).
```julia-repl
julia> plot(scans[1])
julia> plot(chr)
```
"""
module plots

using Plots, RecipesBase   # used for plotting


using MassJ:MScontainer
using MassJ:IonCurrent
using MassJ:maxic
using MassJ:MSscans
using MassJ:YieldCurve
using MassJ:stdev, sem

"""
    normalisation(ms::MassJ.MScontainer)
Normalization function for plotting mass spectra in relative intensity.
"""
function normalisation(ms::MScontainer)
    factor = 100. / ms.basePeakIntensity
    return ms.int .* factor
end

"""
    normalisation(cr::MassJ.IonCurrent)
Normalization function for plotting ion-current traces in relative intensity.
"""
function normalisation(cr::IonCurrent)
    factor = 100. / maxic(cr)
    return cr.ic .* factor
end

"""
    scaling(cr::MassJ.IonCurrent)
Scaling function to display retention times of chromatograms in minutes instead of seconds.
Other axes (drift time, compensation voltage) are returned unscaled.
"""
function scaling(cr::IonCurrent)
    return cr.axis === :rt ? cr.x ./ 60. : cr.x
end


"""
    g(ms::MassJ.MSscans; method = :relative, band = :none)
Allows plotting mass spectra directly. The default relative-intensity plotting
may be changed by setting `method = :absolute`. Set `band = :sem` or `band = :std`
to draw a ±1-σ uncertainty ribbon (standard error of the mean, or sample standard
deviation) on an averaged spectrum; it is ignored when no replicate statistics are
present.
"""
@recipe function g(ms::MSscans; method = :relative, band = :none)
    seriestype --> :path
    seriescolor --> :red
    label --> ""
    xlabel --> "m/z"
    if method == :relative
        factor = 100. / ms.basePeakIntensity
        ylabel --> "Intensity (%)"
    else  # :absolute
        factor = 1.0
        ylabel --> "Intensity (a.u.)"
    end
    if (band === :std || band === :sem) && !isempty(ms.s)
        e = band === :sem ? sem(ms) : stdev(ms)
        ribbon    --> e .* factor
        fillalpha --> 0.15
    end
    ms.mz, ms.int .* factor
end


"""
    h(cr::MassJ.IonCurrent; method = :relative)
Allows plotting ion-current traces directly. The default relative-intensity plotting may be changed by setting method = :absolute.
"""
@recipe function h(cr::IonCurrent; method = :relative)
    seriestype  --> :path
    seriescolor --> :blue
    fillrange   --> 0
    fillalpha   --> 0.3
    label       --> ""
    xlabel      --> (cr.axis === :rt ? "time (mins)" :
                     cr.axis === :drift ? "drift time" : "compensation voltage")
    if method == :relative
        y = normalisation(cr)
        ylabel  --> "Intensity (%)"
    elseif method == :absolute
        y = cr.ic
        ylabel  --> "Intensity (a.u.)"
    end
    scaling(cr), y
end


"""
    k(yc::MassJ.YieldCurve)
Plot a YieldCurve: one line per peak, x-axis from `yc.x`, legend from `yc.labels`,
x-label from `yc.xlabel`. When `yc.yields_err` carries finite uncertainties,
1-σ ribbons are drawn around each line with `fillalpha = 0.15`.
"""
@recipe function k(yc::YieldCurve)
    seriestype --> :path
    xlabel     --> yc.xlabel
    ylabel     --> "yield (a.u.)"
    label      --> reshape(yc.labels, 1, :)
    if any(isfinite, yc.yields_err)
        ribbon    --> replace(yc.yields_err, NaN => 0.0)
        fillalpha --> 0.15
    end
    yc.x, yc.yields
end


end # submodule
