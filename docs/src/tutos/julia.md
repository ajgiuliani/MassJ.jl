# The Julia language
Julia is an open source programming language designed for scientific and technical computing. This section will give a very brief introduction to the Julia language. 

!!! note

    See also the [Julia](https://julialang.org/) language home page, the [Wikibook](https://en.wikibooks.org/wiki/Introducing_Julia) introduction to Julia or [Julia by example](https://juliabyexample.helpmanual.io).
	
## Installation

Julia binaries are available for various platforms and can be downloaded [here](https://julialang.org/downloads).  Plateform specific instructions may be found [here](https://julialang.org/downloads/platform.html).

## Executing Julia code
Julia code can be executed interactively using the REPL as follow:

```@repl
println("Hello")
```

The `println("Hello")` line can be put in a file, such as `hello.jl` and executed as a script, which will produce the same output.
```
$ julia hello.jl
```

You can also put `#!/usr/local/bin/julia` in the first line of the `hello.jl` file, make it executable (` chmod +x hello.jl`) and execute it like any other executable.

Finally, Julia scripts can be executed within Jupyter Notebooks, see the dedicated section [Jupyter notebooks](jupyter.md).


## Types
Julia is a strongly typed language ( see [the julia wikibook](https://en.wikibooks.org/wiki/Introducing_Julia/Types#Type_hierarchy)) In Julia types are organized in a hierarchy with a tree structure. The root of the tree is the `Any` type.  The `Number` type is a direct child of `Any` and possesses two subtypes: `Complex` and `Real`. The `Real` type has four types: `Integer`, `AbstractFloat`, `Irrational` and `Rational`. The exact number hierarchy is as follow:
```
Number
 ↳ Complex
 ↳ Real
    ↳ AbstractFloat 
	   ↳ BigFloat 
	   ↳ Float64 
	   ↳ Float32 
	   ↳ Float16
    ↳ Integer
	   ↳ BigInt
	   ↳ Bool
	   ↳ Signed
	      ↳ Int128 
		  ↳ Int64 
		  ↳ Int32 
		  ↳ Int16 
		  ↳ Int8
	   ↳ Unsigned
	      ↳ UInt128 
		  ↳ Int64 
		  ↳ UInt32 
		  ↳ UInt16 
		  ↳ UInt8
	↳ Irrational
	↳ Rational
```

## Creating vectors and matrices
A vector is created as follow
```@repl
A = [1, 2, 3]                  # vector
A = range(1, 10, step = 2)   # linearly spaced
A = range(1, 10, length = 5) # linearly spaced
A = rand(10)                 # random with 10 elements
j = 2; k = 2; n = 10;
A = j:k:n                    # from j to n with step size k
```
and similarly for matrices
```@repl
A = [1 2; 3 4]               # matrix
A = rand(2, 2)               # random 2x2 matrix
```

Vector and matrices can be manipulated as follow:
```@jldoctest
transpose(A)                 # Return the transpose of A
A[:]                         # Flatten matrix A (convert matrix to vector)
A[2,2]                       # Accessing element at row 2 and colomun 2
A[1:4, :]                    # Accessing specific rows 1 to 4
```

Vector and matrix may be preallocated like this:
```@repl
A = rand(5)                  # a vector / matrix
B = similar(A)               # an emply vector / matrix similar to A
```

## Operations
The following present a few example of operations:
```@jldoctest
julia> dot(A,B)                     # dot product between A and B 
julia> A .* B                       # Element wise multiplication
julia> A * B                        # Matrix multiplcation
julia> norm(A)                      # Euclidian norm
julia> sum(A, dims = 1)             # sum over each column
julia> sum(A, dims = 2)             # sum over each rows
```

## Loops
A loop in Julia can be done like this:
```@jldoctest
julia> for i in 1:N
julia>   #do something
julia> end
```
While loops my be achieved like this:
```@jldoctest
julia> while i <= N
julia>   #do something
julia> end
```
and if / else flow like this:
```@jldoctest
julia> if i < N
julia>   #do something
julia> else 
julia>   #do something else
julia> end
```

## Functions and methods
A function is defined like this:
```@jldoctest
julia> function f(x)
julia>   return x^2
julia> end
```
which can be simplified as:
```@jldoctest
julia> f(x) = x^2
```

Broadcasting a function over a collection or an Array is achieved like this:
```@jldoctest
julia> f(x) = x^2
julia> x = 1:10
julia> f.(x)
```

In Julia, functions that modify their arguments are named with an exclamation mark `!` at the end of their name, such as:
```@jldoctest
julia> function f!(out, x)
julia>   out = x.^2
julia> end
julia> x = rand(10)
julia> y = similar(x)
julia> f!(y, x)
```


## Importing and using Packages
The Julia code is organized into files, modules and packages. A file using the `.jl` extension contains julia code.
Related functions and variable may be gathered in `modules`.  One or more modules may be organized into `packages`. To use a package, it has to be called like this:
```@jldoctest
julia> using A_package
```
If it is not installed, before being used:
```@jldoctest
julia> using Pkg                      # using the Package manager package
julia> Pkg.add("A_package")
julia> using A_package
```
Now every public function which is made available by the package `A_package` is available directly. Private functions have to be called like this: 
```@jldoctest
julia> a_public_function()                   # Calling a public function
julia> A_package.a_private_function()        # Calling private functions
```
The same is true for the variables defined in the packages.



# Jupyter notebooks
[Jupyter](https://jupyter.org/) notebooks are web based documents that may contains both codes, figures and other textual elements (such as equations, links, ...).  Jupyter notebook may be installed easily using Julia:
```
julia> using Pkg
julia> Pkg.add("IJulia")
```
When IJulia is installed, then a notebook may be launch like this:
```
julia> using IJulia
julia> notebook()
```
The `notebook()` function should launch a web browser from which a new notebook  may be started. On each entry of the notebook code, Markdown or text may be inserted. Each line of code may be executed and will eventually return a result.
In the following, the tutorials are organized in [Jupyer notebooks](@ref) form and can be viewed using [nbviewer](http://nbviewer.jupyter.org/).
