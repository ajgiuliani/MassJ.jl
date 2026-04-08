# Properties calculations
From a molecular formula, it is possible to : 
- calculate the molecular masses 
- calculate the isotopic distributions
- simulate a mass spectrum providing and isotopic distribution and a peak width.

## Parsing a chemical formula
The private function [`MassJ.formula`](@ref) takes an input string representing a chemical formula, such as "CH4", counts the number of atoms in the formula and returns a [dictionary](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1). 
```julia
MassJ.formula("CH3Br")
Dict("Br" => 1,"C" => 1,"H" => 3)
Dict{String,Int64} with 3 entries:
  "Br" => 1
  "C"  => 1
  "H"  => 3
```
In the previous example, the dictionary had three elements corresponding to the three atomic species found in the "CH3Br" formula. The "Br" key has the value 1, the "C" key has also 1 and the "H" key has a value equals to 3. The output of this function will be used by all the other functions, which will calculate properties from it.
The  [`MassJ.formula`](@ref) function accepts stoichiometric regular molecular formulas as well as more developed forms. For hexane the following entries "C6H14", "CH3CH2CH2CH2CH2CH3" or "CH3(CH2)4CH3" are equivalents. 
```julia
MassJ.formula("C6H14");
Dict("C" => 6,"H" => 14)
MassJ.formula("CH3CH2CH2CH2CH2CH3");
Dict("C" => 6,"H" => 14)
MassJ.formula("CH3(CH2)4CH3");
Dict("C" => 6,"H" => 14)
```
The indices found after a parenthesis are used to multiply the elements inside the parenthesis. If no indices are found after ")", then the group is not multiply.
The parenthesis are also used to specify an isotope, for example, consider ethane isotopically labelled with carbon 13:
```julia
MassJ.formula("CH3(13C)H3");
Dict("C" => 1,"13C" => 1,"H" => 6)
```
These expressions are equivalent: "CH3(13C)H3", "CH3(13CH3)", "C(13C)H6".
The following isotopes are recognized by MassJ.formula:
- 1H   for ``{}^1H_{1}`` (protium)
- 2H   for ``{}^2H_{2}`` (deuterium)
- D    is equivalent to 2H
- 12C  for ``{}^{12}C_{6}``
- 13C  for ``{}^{13}C_{6}``
- 14N  for ``{}^{14}N_{7}``
- 15N  for ``{}^{15}N_{7}``
- 16O  for ``{}^{16}O_{8}``
- 18O  for ``{}^{18}O_{8}``
- 32S  for ``{}^{32}S_{16}``
- 34S  for ``{}^{34}S_{16}``

The `formula` function is also exported as a public function and can be called directly:
```julia
julia> f = formula("C6H14")
Dict("C" => 6,"H" => 14)
```

The elements are stored in a constant dictionary called `MassJ.Elements`. Each key of the `MassJ.Elements` points to an `Array` of [`MassJ.Isotope`](@ref), which is a structure used to stored the different properties of the isotopes:
```julia
struct Isotope
    m::Float64           # mass
    f::Float64           # natural frequency
    logf::Float64        # logarythm of the natural frequency
    Z::Int               # atomic number
    A::Int               # mass number
    active::Bool         # is radioactive
end
```
The isotopes are sorted by natural frequency. Hence, for instance, the first sulfur isotope is ``{}^{32}S_{16}`` with a natural frenquency of 0.995 followed by  ``{}^{34}S_{16}`` with 0.043, etc.
```julia
julia> E = MassJ.Elements["S"]
4-element Array{MassJ.Isotope,1}:
 MassJ.Isotope(31.9720711741, 0.9498500119990401, -0.051451188958515866, 16, 32, false)
 MassJ.Isotope(33.96786703, 0.04252059835213182, -3.157766653355949, 16, 34, false)
 MassJ.Isotope(32.9714589101, 0.00751939844812415, -4.890269137820559, 16, 33, false)
 MassJ.Isotope(35.9670812, 0.00010999120070394368, -9.115110188972029, 16, 36, false)
```
The first element of the array, is the most naturally abundant isotope. The properties of the isotopes may be accessed by:
```julia
julia> E = MassJ.Elements["S"];
julia> E[1]                            # returns the most abundant isotope of sulfur
MassJ.Isotope(31.9720711741, 0.9498500119990401, -0.051451188958515866, 16, 32, false)
julia> E[2].f                          # returns natural abundance of the second most abundant istotope of sulfur.
0.04252059835213182
```


## Molecular masses
The public function [`masses`](@ref), return a dictionary with three keys "Monoisotopic", "Average" and "Nominal", defined as:

| Mass        | Element                                          | Molecule or ion                                                |
|-------------|--------------------------------------------------|----------------------------------------------------------------|
| Nominal     | Mass number of the most abundant stable isotope  | Sum of the nominal masses of the most abundant stable isotopes |
| Monoistopic | Atomic mass of the most abundant stable isotope  | Sum of the atomic masses of the most abundant stable isotopes  |
| Average     | Mass average function of the relative abundances | Sum of the average atomic weights of the constituents          |

The masses for hexane are calculated by:
```julia
julia> m = masses("CH3(CH2)4CH3")
Dict("C" => 6,"H" => 14)
Dict{String,Float64} with 3 entries:
  "Monoisotopic" => 86.1096
  "Average"      => 86.178
  "Nominal"      => 86.0
```
And the specific masses may be accessed by:
```julia
julia> m["Average"]
86.178
julia> m["Monoisotopic"]
86.10955045178
julia> m["Nominal"]
86.0
```

## Isotopic distributions
The public function [`isotopic_distribution`](@ref) calculates the isotopic distribution of the given formula. The function takes the following arguments:
- formula: `String`
- target probability: `Real` number 
- charge: optional argument `Int`, by default = +1
- tau: optional `Real` number set by default to 0.1
- the elements dictionary: set by default to `MassJ.Elements`.

The calculations is based on the implementation of the [`isospec`](https://doi.org/10.1021/acs.analchem.6b01459) algorithm.  Briefly, the algorithm search for the small set of [`isotopologues`](https://en.wikipedia.org/wiki/Isotopologue) for which the total abundance is equal to the target probability. The calculation return a vector with all the configurations found. The first column gives the masses, the second column the probabilities of the configurations, and the following columns gives the configurations, such as:
```julia
julia> I = isotopic_distribution("C6H14", 0.9999, charge = +1)
Dict("C" => 6,"H" => 14)
5×6 Array{Union{Float64, Int64, String},2}:
   "Masses"   "Probability"   "12C"   "13C"    "1H"   "2H"
 86.1096     0.935476        6       0       14      0
 87.1129     0.0612122       5       1       14      0 
 88.1163     0.00166891      4       2       14      0 
 87.1158     0.00151559      6       0       13      1
```
 

The resulting array is intend to be readable easily, and thus the first line contains the descriptors or the columns. It may also be easily exported to a delimited file.

## Simulated mass spectra
The result of an isotopic distribution calculation may be convoluted with a peak shape to produce a simulated mass spectrum. Such an operatio, is achieved by the [`simulate`](@ref) function.
The function takes the following arguments:
- the `Array`resulting from the `isotopic_distribution` calculation
- the width in Dalton of the peak shape
- the model used for the peak shape. By default it is set to `model=:gauss`, but may be `model=:lorentz` or `model=:voigt`.
- An `Int` value for the number of points in the mass spectrum. By default set to `Npoints=1000`.

Hence, the simulated istopic distribution of haxane can be obtained and plotted like this:
```julia
julia> I = isotopic_distribution("C6H14", 0.9999, charge = +1);
Dict("C" => 6,"H" => 14)
julia> sim = simulate(I, 0.05, model= :lorentz);
julia> Using Plots
julia> plot(sim)
```


See [Plotting](@ref) for more information on how to plot mass spectra.

