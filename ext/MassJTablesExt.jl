module MassJTablesExt

using MassJ: YieldCurve
import Tables

# Wide column table for a YieldCurve: one row per x value (file), columns
#   x, <peak label>, <peak label>_err, …, TIC, TIC_err, file
# A NamedTuple of vectors is itself a valid Tables column table, so this is all
# the Tables interface needs.
function _columntable(yc::YieldCurve)
    names = Symbol[:x]
    cols  = Any[yc.x]
    @inbounds for (j, lbl) in enumerate(yc.labels)
        push!(names, Symbol(lbl));         push!(cols, yc.yields[:, j])
        push!(names, Symbol(lbl, :_err));  push!(cols, yc.yields_err[:, j])
    end
    push!(names, :TIC);     push!(cols, yc.tic)
    push!(names, :TIC_err); push!(cols, yc.tic_err)
    push!(names, :file);    push!(cols, yc.files)
    return NamedTuple{Tuple(names)}(Tuple(cols))
end

Tables.istable(::Type{<:YieldCurve})      = true
Tables.columnaccess(::Type{<:YieldCurve}) = true
Tables.columns(yc::YieldCurve)            = _columntable(yc)
Tables.rowaccess(::Type{<:YieldCurve})    = true
Tables.rows(yc::YieldCurve)               = Tables.rows(_columntable(yc))

end # module
