"""
Main module for `MassJ.jl`-- A Julia package to load and process mass spectrometry data.

"""
module MassJ


using Statistics           # used for Pearson's correlation calculation, mean
using LsqFit               # used for curve fitting
using FFTViews             # used for fft
using Interpolations: Interpolations, LinearInterpolation, Line  # used for interpolation
using LinearAlgebra        # used for matrix operation 
using LightXML, Codecs     # used for mzXML file import
using Unicode              # used for file io
using Combinatorics        # used for factorial & permutations
using DataStructures       # used for PriorityQueue
using Printf               # used for @sprintf
using Libz                 # used for zlib compression (mzxml import)
using ProgressMeter        # used in deconvolution
using DelimitedFiles       # used in txt file import
using Parameters           # used to construct Charge and Mass struct

import Base: +, -, *, /


include("types.jl")
include("Io.jl")
include("msscans.jl")
include("mzxml.jl")
include("mzml.jl")
include("mgf.jl")
include("msp.jl")
include("imzml.jl")
include("process.jl")
include("extract.jl")
include("utilities.jl")
include("predicates.jl")
include("isotopes.jl")
include("deconvolution.jl")
include("txt.jl")


# Submodules
# ----------

include("plots.jl")



end # module
