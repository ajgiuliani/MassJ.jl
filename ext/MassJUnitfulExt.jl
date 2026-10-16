module MassJUnitfulExt

using MassJ
using MassJ: MSscans, IonCurrent
using Unitful: @u_str

# Physical-axis quantities of a scan, tagged with units. m/z (a dimensionless
# ratio) and intensity (arbitrary units) are intentionally left untagged. Drift
# time is omitted because the field may hold a drift time (ms) or a reduced
# mobility 1/K₀ depending on the instrument — a single unit would be wrong.
function MassJ.withunits(scan::MSscans)
    return (retention_time       = scan.rt .* u"s",
            collision_energy      = scan.collisionEnergy .* u"eV",
            compensation_voltage  = scan.compensationVoltage .* u"V")
end

# An ion-current trace's abscissa, tagged according to its axis (intensity left
# untagged). The :drift axis is returned unitless for the reason above.
function MassJ.withunits(ic::IonCurrent)
    x = ic.axis === :rt ? ic.x .* u"s" :
        ic.axis === :cv ? ic.x .* u"V" :
        ic.x
    return (x = x, ic = ic.ic)
end

end # module
