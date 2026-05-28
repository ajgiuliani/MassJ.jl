

"""
    filter(scans::Vector{MSscans}, argument::Level{<:Int})
Search for scans matching the argument MS level and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Level{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.level[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Level{<:AbstractVector})
Search for scans matching the argument MS levels and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Level{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.level[1] == i
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Precursor{<:Real})
Search for scans matching the argument precursor mz and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Precursor{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.precursor[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::Precursor{<:AbstractVector})
Search for scans matching the argument precursors mz and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Precursor{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.precursor[1] == i
                push!(subindex, elem.num[1])
            end
        end        
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Activation_Energy{<:Real})
Search for scans matching the argument activation energy and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Activation_Energy{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.collisionEnergy[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Activation_Energy{<:AbstractVector})
Search for scans matching the argument activation energies and returns a list of the corresponding indexes

"""
function filter(scans::Vector{MSscans}, argument::Activation_Energy{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for elem in scans
            if elem.collisionEnergy[1] == i
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Activation_Method{<:String})
Search for scans matching the argument activation method and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Activation_Method{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.activationMethod[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Activation_Method{<:AbstractVector})
Search for scans matching the argument activation methods and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Activation_Method{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.activationMethod[1] == i
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{MSscans}, argument::Polarity{<:AbstractVector})
Search for scans matching the argument polarities and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Polarity{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.polarity[1] == i
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::Polarity{<:String})
Search for scans matching the argument polarity and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Polarity{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.polarity[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::Scan{<:Int})
Search for scans matching the argument scan num and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Scan{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.num[1] == argument.arg
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::Scan{<:AbstractVector})
Search for scans matching the argument scan nums and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::Scan{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for elem in scans
            if elem.num[1] == i
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::RT{<:Real}) 
Search for scans matching the argument retention time and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::RT{<:Real}) 
    subindex = Set{Int}()
    rt = retention_time(scans)
    index = num2pnt(rt, argument.arg)
    for elem in scans
        if elem.num[1] == index
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::RT{<:AbstractVector})
Search for scans matching the argument retention time in the specified range and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::RT{<:AbstractVector})
    subindex = Set{Int}()
    rt = retention_time(scans)
    for i =1:2:length(argument.arg)
        index_low  = num2pnt(rt, argument.arg[i])
        index_high = num2pnt(rt, argument.arg[i+1])
        for elem in scans
            if  index_low <= elem.num[1] <= index_high
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::RT{<:AbstractVector{<:AbstractVector} } )
Search for scans matching the argument retention time in the specified ranges and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::RT{<:AbstractVector{<:AbstractVector} } )
    subindex = Set{Int}()
    rt = retention_time(scans)
    for el  in argument.arg
        index_low  = num2pnt(rt, el[1])
        index_high = num2pnt(rt, el[2])
        for elem in scans
            if  index_low <= elem.num[1] <= index_high
                push!(subindex, elem.num[1])
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{MSscans}, argument::IC{<:AbstractVector})
Search for scans matching the argument total ion current within the specified ranges and returns a list of the corresponding indexes
"""
function filter(scans::Vector{MSscans}, argument::IC{<:AbstractVector})
    subindex = Set{Int}()
    for elem in scans
        if  argument.arg[1] <= elem.tic <= argument.arg[2]
            push!(subindex, elem.num[1])
        end
    end
    return subindex
end



### Extraction of the ion current according to the selected filters and method

"""
    extracted_chromatogram(scans::Vector{MSscans}, indices::Vector{Int},method::MethodType)
Returns the extracted chromatogram for input Array of MSscan according to the selected method and for set of scan num as input
"""
function extracted_chromatogram(scans::Vector{MSscans}, indices::Vector{Int},method::MethodType)
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    if method isa BasePeak
        for i = 1:length(indices)
            push!(xrt, scans[indices[i]].rt[1])
            push!(xic, scans[indices[i]].basePeakIntensity)
        end
    elseif method isa ∆MZ
        mz1 = convert(Float64, method.arg[1] - method.arg[2] )  # mz - ∆mz
        if(mz1 < 0.0)
            return ErrorException("Bad mz ± ∆mz values.")
            #error("Bad mz ± ∆mz values")
        end
        mz2 = convert(Float64, method.arg[1] + method.arg[2] ) # mz + ∆mz
        for i = 1:length(indices)
            value = add_ion_current(scans[indices[i]].mz, scans[indices[i]].int, mz1, mz2)
            push!(xrt, scans[indices[i]].rt[1])
            push!(xic, value)
        end

    elseif method isa MZ
        mz1 = convert(Float64,method.arg[1])
        mz2 = convert(Float64,method.arg[2])
        for i = 1:length(indices)
            value = add_ion_current(scans[indices[i]].mz, scans[indices[i]].int, mz1, mz2)
            push!(xrt, scans[indices[i]].rt[1])
            push!(xic, value)
        end
    else
        for i = 1:length(indices)
            push!(xrt, scans[indices[i]].rt[1])
            push!(xic, scans[indices[i]].tic)
        end

    end
    return IonCurrent(xrt, xic; axis = :rt)
end




### Calculation of the composite MSscans composite mass spectrum according to the selected filters and method


"""
    composite_spectra(scans::Vector{MSscans}, indices::Vector{Int}, stats::Bool)
Returns the average MSscans for input Array of MSscan and according to the input scan num. Calculation of variance is controlled by the stats Boolean variable.
"""
function composite_spectra(scans::Vector{MSscans}, indices::Vector{Int}, stats::Bool)
    if stats == false
        result = scans[indices[1]]
        for i = 2:length(indices)
            result += scans[indices[i]]
        end
        return result / length(indices)
    elseif stats == true
        result = scans[indices[1]]
        for i = 2:length(indices)
            result = avg(result, scans[indices[i]])
        end
        return standard_deviation(result, length(indices))
    end
end

