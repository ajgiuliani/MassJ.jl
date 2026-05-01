"""
Isotopic distributions.

Independent reimplementation of the IsoSpec algorithm. Reference:
  Łącki, M. K.; Startek, M.; Valkenborg, D.; Gambin, A.
  "IsoSpec: Hyperfast Fine Structure Calculator." Anal. Chem. 2017, 89, 3272.
"""

include("isotopes-data.jl")

# User Interface.
# ---------------


export isotopic_distribution, masses, simulate, formula


"""
    masses(input::String)
Calculates the average, monoisotopic and nominal masses for the chemical formula given as an input. The result is returned in a dictionary with the following entries: "Monoisotopic", "Average" and "Nominal".
# Examples
```julia-repl
julia> masses("C254 H377 N65 O75 S6")
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
Dict{String,Float64} with 3 entries:
  "Monoisotopic" => 5729.6
  "Average"      => 5733.55
  "Nominal"      => 5727.0
```
"""
function masses(input::String)
    f = formula(input)
    masses(f)
end

"""
    masses(f::Dict{String,Int})
Calculates the average, monoisotopic and nominal masses for the chemical formula dictionary, such as produced by MassJ.formula. The result is returned in a dictionary with the following entries: "Monoisotopic", "Average" and "Nominal".
# Examples
```julia-repl
julia> masses("C254 H377 N65 O75 S6")
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
Dict{String,Float64} with 3 entries:
  "Monoisotopic" => 5729.6
  "Average"      => 5733.55
  "Nominal"      => 5727.0
```
"""
function masses(f::Dict{String,Int})
    m = Dict{String, Float64}()

    # exact mass
    avg = 0.0
    mono = 0.0
    nomi = 0.0
    for (key,val) in f
        avg += AverageMasses[key] * val
        mono += first(Elements[key]).m * val
        nomi += Int( map(round,first(Elements[key]).m) ) * val
    end
    m["Average"] = avg
    m["Monoisotopic"] = mono
    m["Nominal"] = nomi
    m
end


"""
    simulate(I::Array{Union{Float64, Int, String}}, ∆mz::Real; model::Symbol=:gauss, Npoints::Int=1000)
From an isotopic distribution and a peak width returns a mass spectrum (MSScan). The number of points of the resulting mass spectrum is passed as an optional argument. Peak shape are :gauss (default), :lorentz, :voigt.
# Examples
```julia-repl
julia>  a = simulate(I, 0.4)
MassJ.MSscan(1, 0.0, 30898.192348114364, [5727.102517458742 ..., "", "", 0.0)
```
"""
function simulate(I::Array{Union{Float64, Int, String}}, ∆mz::Real; model::Symbol=:gauss , Npoints::Int=1000)
    r =  maximum(I[2:end,1]) - minimum(I[2:end,1])
    start =  minimum(I[2:end,1]) - r/2.
    stop  =  maximum(I[2:end,1]) + r/2.
    r = stop - start
    
    p  = zeros(Npoints)
    mz = zeros(Npoints)
    peak = zeros(Npoints)  
    
    mz = [start + (x * r / length(mz)) for x=1:length(mz)]
    
    for (i,val) in enumerate(I[2:end, 1])
        ind = num2pnt(mz, val)
        x0 = mz[ind]
        y0 = I[2:end,2][i]
        param = [∆mz, x0, y0, 0.0]
        if model == :gauss
            peak = [gauss(x,param) for x in mz]
        elseif model == :lorentz
            peak = [lorentz(x,param) for x in mz]
        elseif model == :voigt
            peak = [voigt(x,param) for x in mz]
        else
            ErrorException("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
        end            
        p = p .+ peak
    end
    p /= maximum(p)
    p *= 100.0
    basePeakMz = mz[num2pnt(p, 100.0)]
    MSscan(1, 0.0, sum(p), mz, p, 0, basePeakMz, 100.0, 0.0, "", "", 0.0)
end

"""
    isotopic_distribution(input::String, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{MassJ.Isotope,1}} = MassJ.Elements)
Calculates the isotopic distribution of input formula for which the overall probabilities equals p_target using the isospec algorithm. The charge state is entered as an optional argument. The peaks detection threshold tau is by default set to 10%.
# Examples
```julia-repl
julia> isotopic_distribution("C254 H377 N65 O75 S6", 0.5)
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
9×15 Array{Union{Float64, Int, String},2}:
     "Masses"   "Probability"     "12C"   "13C"   "32S"   "34S"   "33S"   "36S"    "14N"   "15N"    "16O"   "18O"   "17O"     "1H"   "2H"
 5731.61       0.112302        252       2       6       0       0       0       65       0       75       0       0       377      0    
 5732.61       0.102878        251       3       6       0       0       0       65       0       75       0       0       377      0    
 5730.6        0.0814037       253       1       6       0       0       0       65       0       75       0       0       377      0    
 5733.61       0.0704028       250       4       6       0       0       0       65       0       75       0       0       377      0    
 5734.62       0.0383896       249       5       6       0       0       0       65       0       75       0       0       377      0    
 5733.6        0.0301637       252       2       5       1       0       0       65       0       75       0       0       377      0    
 5729.6        0.0293871       254       0       6       0       0       0       65       0       75       0       0       377      0    
 5734.61       0.0276323       251       3       5       1       0       0       65       0       75       0       0       377      0    
``julia> isotopic_distribution("C254 H377 N65 O75 S6", 0.5, charge = +7)
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
9×15 Array{Union{Float64, Int, String},2}:
    "Masses"   "Probability"     "12C"   "13C"   "32S"   "34S"   "33S"   "36S"    "14N"   "15N"    "16O"   "18O"   "17O"     "1H"   "2H"
 818.801      0.112302        252       2       6       0       0       0       65       0       75       0       0       377      0    
 818.944      0.102878        251       3       6       0       0       0       65       0       75       0       0       377      0    
 818.658      0.0814037       253       1       6       0       0       0       65       0       75       0       0       377      0    
 819.088      0.0704028       250       4       6       0       0       0       65       0       75       0       0       377      0    
 819.231      0.0383896       249       5       6       0       0       0       65       0       75       0       0       377      0    
 819.086      0.0301637       252       2       5       1       0       0       65       0       75       0       0       377      0    
 818.514      0.0293871       254       0       6       0       0       0       65       0       75       0       0       377      0    
 819.23       0.0276323       251       3       5       1       0       0       65       0       75       0       0       377      0    
```
"""
function isotopic_distribution(input::String, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{MassJ.Isotope,1}} = MassJ.Elements)
    f = formula(input)
    println(f)
    isotopic_distribution(formula(input), p_target, charge = charge, tau = tau)
end

function isotopic_distribution(form::Dict{String,Int}, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{MassJ.Isotope,1}} = MassJ.Elements)
    c_form = collect(form)
    sort_formula(x) = 
        begin 
            if length(Elements[x]) > 1
                return Elements[x][2].f * form[x]
            else
                return 0.0
            end
        end
    sort!(c_form, by = x -> sort_formula(x[1]), rev=true)    
    I = isospec(form, c_form, p_target, tau, Elements)    
    cI = collect(I)
    
    N = 0
    for (key, val) in c_form
        N += length(Elements[key])
    end
    N += 2                              # columns
    M = length(I)                       # rows
    result = Matrix{Union{String, Float64, Int}}(undef, M+1,N)
    labels = Vector{String}(undef, N)
    labels[1] = "Masses"
    labels[2] = "Probability"
    for i=1:M
        temp_dist = Dict{String,Array{Int,1}}()
        j = 1
        for k = 1:length(c_form)
            elem = c_form[k][1]
            temp_vect = [ cI[i][1][j:j+length(Elements[elem])-1]... ]
            result[i,j+2:j+2+length(Elements[elem])-1] = temp_vect
            j = j + length(Elements[elem])
            temp_dist[elem] = temp_vect
        end
        masses =  isotopologue_mass(c_form, temp_dist, Elements)
        result[i,1] = masses / abs(charge)
        result[i,2] = cI[i][2]
    end
    j = 3
    for i=1:length(c_form)
        elem = c_form[i][1]
        isot = length(Elements[elem])
        for m=1:length(Elements[elem])
            labels[j] = string(Elements[elem][m].A) * elem
            j+=1
        end
    end
    result[M+1,:] = labels
    result[end:-1:1,:]
end


"""
    formula(formula::String)
Private function that reads the input chemical formula and sorts the atoms. It returns a dictionary in which the different entries represent the atoms and the values are the number of times the atoms have been found in the formula.
"""
function formula(formula::String)
    if !occursin(r"[A-Z-(]", formula[1:1])                                        # is Formula starting with capital ?
        error("Invalid chemical formula")
    end    
    form_dict = Dict{String,Int}() 
    i = 1
    j = 1
    while i <= length(formula)
        if occursin(r"[A-Z]", formula[i:i])                                     # is formula[i] is a capital
            temp = formula[i:i]
            prev = (i:i)
            if i < length(formula) 
                if occursin(r"[a-z]", formula[i+1:i+1])                         # is the capital is followed by a lower-case ? 
                    temp = formula[i:i+1]                 
                    prev = (i:i+1)
                end
            end
            if haskey(MassJ.Elements, temp)                                   # valid element --> looking for numbers 
                if last(prev) == length(formula)                                # last element
                    if !haskey(form_dict, formula[prev])                           # already found ?
                        form_dict[formula[prev]] = 1                               # no --> indice = 1
                    else
                        form_dict[formula[prev]] = form_dict[formula[prev]] + 1 # yes --> add 1 to indice
                    end 
                elseif last(prev) < length(formula)                             # not last element
                    if occursin(r"[0-9]", formula[last(prev):last(prev)+1])        # indice after element ?               
                        next = findnext(r"\d+", formula, last(prev)+1)
                        if !haskey(form_dict, formula[prev])                          # already found ?
                            form_dict[formula[prev]] = parse(Int,formula[next])  
                        else 
                            form_dict[formula[prev]] = form_dict[formula[prev]] + parse(Int,formula[next])
                        end
                    else                                                       # no indice after element
                        if !haskey(form_dict, formula[prev])                          # already found ?
                            form_dict[formula[prev]] = 1                                # no --> indice = 1
                        else
                            form_dict[formula[prev]] = form_dict[formula[prev]] + 1 
                        end    
                    end    
                end
            else
                error("Invalid chemical element: $temp")
            end
        elseif occursin(r"[(]", formula[i:i])                                  # formula contains a ( --> isotopes or group
            j = i+1
            mult = 1
            l = j
            while !occursin(r"[)]", formula[l+1:l+1])         # find closing ) and check for an indice --> mult factor
                l+=1
            end
            if l+2 > length(formula)                       # if ) is the last char of the formula, mult = 1
                mult = 1
            else
                if occursin(r"[1-9]", formula[l+2:l+2])   #  indice found after )
                    k = l+2
                    if k < length(formula)
                        while occursin(r"[0-9]", formula[k:k])
                            k += 1
                        end
                        mult = parse(Int,formula[l+2:k-1])
                    else
                        mult = parse(Int,formula[l+2:k])
                    end
                else
                    mult = 1
                end
            end
            #println("mult :$mult")
            while !occursin(r"[)]", formula[j:j])       # repeat until the end
                l = j
                if occursin(r"[0-9]", formula[j:j])                      # number found after ( --> isotope
                    while occursin(r"[0-9]", formula[l:l])
                        l += 1
                    end                    
                    prev = formula[j:l]
                    if haskey(MassJ.Elements, prev)               # element found
                        j=l
                        if occursin(r"[)]", formula[j+1:j+1])        # no indice after element --> indice = 1
                            next = "1"
                        elseif occursin(r"[1-9]", formula[j+1:j+1])
                            l = j+1
                            while occursin(r"[0-9]", formula[l:l])
                                l+=1
                            end
                            next = formula[j+1:l-1]
                            j = l-1
                        elseif occursin(r"[A-Z]", formula[j+1:j+1])    # found following element --> storing previous one
                            next = "1"
                        end
                        if !haskey(form_dict, prev)
                            form_dict[prev] = parse(Int,next) * mult                                # no --> indice = 1
                        else
                            form_dict[prev] = form_dict[prev] +  (parse(Int,next) * mult) 
                        end                                
                    else
                        error("Invalid chemical element: $prev")
                    end
                elseif occursin(r"[A-Z]", formula[j:j])               # letter found after ( --> chemical group
                    l = j               
                    prev = formula[j:j]
                    next = ""
                    if occursin(r"[a-z]", formula[j+1:j+1])                 # is the capital is followed by a lower-case ? 
                        prev = prev * formula[j+1:j+1]
                        j += 1
                        println(prev)
                        if occursin(r"[)]", formula[j+1:j+1])               # capital letter followed by ) indice -> 1
                            next = "1"
                        elseif occursin(r"[1-9]", formula[j+1:j+1])            # finding indices
                            l = j+1
                            while occursin(r"[0-9]", formula[l:l])
                                l += 1
                            end
                            next = formula[j+1:l-1]
                            j = l-1
                        elseif  occursin(r"[A-Z]", formula[j+1:j+1])        # is followed by another element
                            next = "1"
                        end
                    elseif occursin(r"[1-9]", formula[j+1:j+1])            # finding indices
                        l = j+1
                        while occursin(r"[0-9]", formula[l:l])
                            l += 1
                        end
                        next = formula[j+1:l-1]
                        j = l-1
                    elseif occursin(r"[A-Z]", formula[j+1:j+1])            # found another element next
                        next = "1"
                    elseif occursin(r"[)]", formula[j+1:j+1])               # last element before )
                        next = "1"
                    end
                    if haskey(MassJ.Elements, prev)                             # element found

                        if !haskey(form_dict, prev)                           # already found ?
                            form_dict[prev] = parse(Int,next) * mult         # no --> indice
                        else
                            form_dict[prev] = form_dict[prev] +  (parse(Int,next) * mult) # yes --> add 1 to indice
                        end
                    else
                        error("Invalid chemical element: $prev")
                    end
                end
                j+=1
            end
            i = j
        end
        i += 1
    end

    form_dict
end



function isotopologue_probability(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{MassJ.Isotope,1}})
    p = 1.0                                                                                        # isotopologue probability
    Natoms = sum(collect(values(formula)))
    if Natoms < 20
        return low_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements)
    else

        return high_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements)
    end
    
end

function low_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{MassJ.Isotope,1}})
    p = 1.0
    for (key, val) in formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        
        for i=1:length(distri[key])
            p = p * E[i].f^distri[key][i]
        end
        # getting coefficients
        coefs = factorial(BigInt(val)) 
        for el in distri[key]
            coefs /= factorial(BigInt(el))
        end
        p *= coefs
    end
    p = Float64(p)                            # returns mass and probability for the input isotopologue
end

function high_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{MassJ.Isotope,1}})
    # log (factorial(val) ) - Somme log(factorial(distri[key])) + somme(distri[key] * log(E[i].f) )
    log_p = 0.0
    for (key,val) in formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        
        for i=1:length(distri[key])
            log_p = log_p + distri[key][i]*E[i].logf
        end
        if val <= 256
            log_p = log_p + log_factorial[val]
        else
            log_p  = log_p + stirling(val)
        end
        for i in eachindex(distri[key])
            el = distri[key][i]
            if el != 0
                if el <= 256
                    log_p = log_p - log_factorial[el]
                else
                    log_p = log_p - stirling(el)
                end
            end
            
        end
    end
    exp(log_p)
end

function stirling(n::Int)
    x = n + 1
    return (x - 0.5)*log(x) - x + 0.5*log(2*pi) + 1.0/(12.0*x)
end


function isotopologue_mass(c_formula::Array{Pair{String,Int},1}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{MassJ.Isotope,1}})
    m = 0.0                                      # isotopologue mass
    for (key, val) in c_formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        for i=1:length(distri[key])
            m = m + E[i].m * distri[key][i]
        end
    end
    m                                      # returns the mass of the input isotopologue
end


function most_probable_isotopologue(formula::Dict{String,Int}, Elements::Dict{String,Array{MassJ.Isotope,1}})
    distri=Dict{String,Array{Int,1}}()        
    # enter hill climbing
    for (key, val) in formula
        E = Elements[key]
        if length(E) > 1
            # # start with the mean of the multimodal distribution
            P0 = Int.(map(round, val .* [el.f for el in E]))
            counts = 1
            while sum(P0) != val               # checking for mismatch between distri & formula
                for i=length(P0):-1:1
                    if P0[i] > 0               # found a non negative element                 
                        if sum(P0) > val       
                            P0[i] -= 1
                        elseif sum(P0) < val   
                            P0[i] +=1
                        end
                    end
                end
                if counts < 10
                    break
                end
                counts += 1
            end
            # Objective function isotopologue_probability(Dict(key => val), x::AbstractArray)
            f(x) = isotopologue_probability(Dict(key => val), Dict(key => x), Elements)
            
            # call hill_climbing(P, f) --> returns P
            P = hill_climbing(P0, f)
            distri[key] = P
        else
            # element with single isotopes
            distri[key] = Int.(val .* [el.f for el in E])
        end
    end
    distri
end


function hill_climbing(P::Vector{Int}, f)
    max_count = 5000
    counts = 0
    Pabs = copy(P)
    fPabs = f(Pabs)
    # find all neighbours
    n = fill(0, length(P))
    if length(P) <= 1
        return Pabs
    end
    n[1] = -1
    n[2] = +1
    neighbours = unique(permutations(n))
    temp = similar(P)
    while counts < max_count
        counts += 1
        imax = 0
        fmax = fPabs
        for i in eachindex(neighbours)   # find neighbours that increase the function f
            temp .= P .+ neighbours[i]
            if !any(<(0), temp)                 # no negative values
                ft = f(temp)
                if ft >= fmax
                    fmax = ft
                    imax = i
                end
            end
        end
        if imax != 0
            P = P .+ neighbours[imax]
            fP = fmax
        else
            break
        end

        if fP > fPabs
            Pabs = copy(P)
            fPabs = fP
        else
            break
        end

    end
    if counts >= max_count
        @warn("Couldn't find the most probable subisotopologue.")
    end
    return Pabs
end


function subgenerator!(formula::Dict{String,Int}, D_max::Dict{String,Array{Int,1}},  V::PriorityQueue{Vector{Int},Float64}, τ::Real, Elements::Dict{String,Array{MassJ.Isotope,1}})
    key = first(formula)[1]
    val = first(formula)[2]
    E = Elements[key]
    P_max = D_max[key]
    fP_max = isotopologue_probability(formula, D_max, Elements)

    elem_form = Dict(key => val)
    f(x) = isotopologue_probability(elem_form, Dict(key => x), Elements)

    if isempty(V)
        V[P_max] = f(P_max)
    end
    n = fill(0, length(P_max))
    if length(P_max) > 1                                       # ensure there at least 2 isotopes for this element
        n[1] = -1
        n[2] = +1
        neighbours = unique( permutations(n) )
        threshold = τ * fP_max
        # BFS from existing entries in V
        PQ = PriorityQueue{Vector{Int},Float64}()
        for (k, v) in V
            PQ[k] = v
        end
        while !isempty(PQ)
            beta = dequeue!(PQ)
            for i in eachindex(neighbours)
                temp = beta .+ neighbours[i]
                if !any(<(0), temp)
                    if !haskey(V, temp)                      # if not in V
                        ftemp = f(temp)
                        if ftemp >= threshold              # if prob above threshold
                            V[temp] = ftemp                 # then keep configuration
                            PQ[temp] = ftemp
                        end
                    end
                end
            end
        end
    end
    V
end

function ordered_isotopologue(c_formula::Array{Pair{String,Int},1}, α::Dict{String,Array{Int,1}}, V::Vector{PriorityQueue{Vector{Int},Float64}}, I::PriorityQueue{Vector{Int},Float64}, tau::Real, Elements::Dict{String,Array{MassJ.Isotope,1}} )
    S = Vector{PriorityQueue{Vector{Int},Float64}}(undef,0)

    i = 1
    for (key, val) in c_formula
        form = Dict(key => c_formula[i][2])
        temp = Dict(key => α[key])
        subgenerator!(form, temp, V[i], tau, Elements)
        push!(S, V[i])
        i += 1
    end
    I = S[1]
    for i=2:length(S)
        cI = collect(I)
        I = PriorityQueue{Vector{Int},Float64}()
        cSi = collect(S[i])
        for el in cI
            for j=1:length(cSi)
                conf = vcat(el[1], cSi[j][1])
                prob = el[2] * cSi[j][2]
                I[conf] = prob
            end
        end
    end
    V, I
end

function trim!( I::PriorityQueue{Vector{Int},Float64}, p_target::Real, summ::Real)
    # Remove lowest-probability configurations until total <= p_target
    while !isempty(I)
        lowest_val = peek(I)[2]
        remaining = summ - lowest_val
        if remaining <= p_target
            break
        end
        dequeue!(I)
        summ = remaining
    end
    I
end


function isospec(formula::Dict{String,Int}, c_formula::Array{Pair{String,Int},1}, p_target::Real, tau::Real, Elements::Dict{String,Array{MassJ.Isotope,1}} )
    α = most_probable_isotopologue( formula, Elements )
    I = PriorityQueue{Vector{Int},Float64}()
    N_elem = length(formula)
    V = Vector{PriorityQueue{Vector{Int},Float64}}(undef,0)
    for i =1:N_elem
        push!(V, PriorityQueue{Vector{Int},Float64}())
    end
    count = 1
    summ = 0.0
    while true
        s = @sprintf "%i \t tau: %e" count tau                             # used for debugging
        #print(s)                                                          # used for debugging
        
        V,I = ordered_isotopologue(c_formula, α, V, I, tau, Elements)

        summ = 0.0
        for (key,val) in I
            summ += val
        end
        u = @sprintf "\tprob.: %f \t diff.: %e" summ (p_target-summ)        # used for debugging
        #println(u)                                                         # used for debugging

        a = (summ - 1.0)/tau
        tau = (p_target-1.0) / a
        
        count += 1
        if count >= 5
            break
        end
        if summ >= p_target
            break
        end
    end
    trim!(I, p_target, summ)
    I
end
