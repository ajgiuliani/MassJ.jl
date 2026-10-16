using Documenter, MassJ

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://ajgiuliani.github.io/MassJ.jl/stable/",
        assets = ["assets/favicon.ico"],
        # reference.md is the full single-page API (~190 KiB); raise the per-page
        # size limits (bytes) above it, with headroom as the API grows.
        size_threshold_warn = 300_000,
        size_threshold = 500_000,
    ),
    modules = [MassJ],
    sitename = "MassJ.jl",
    authors = "Alexandre Giuliani",
    warnonly = true,
    checkdocs = :exports,
    pages = [
        "Home"                             => "index.md",
        "Tutorials"                        => Any[
            "The Julia language"           => "tutos/julia.md",
            "Jupyer notebooks"             => "tutos/jupyter.md",
            "MassJ"                        => "tutos/MassJ.md",
        ],
        "Manual"                           => Any[
            "Introduction"                 => "man/introduction.md",
            "Public elements"              => "man/public.md",
            "Data types"                   => "man/types.md",
            "File Information"             => "man/information.md",
            "Importing data"               => "man/importing.md",
            "Exporting data"               => "man/exporting.md",
            "Combining and filtering data" => "man/filtering.md",
            "Processing"                   => "man/processing.md",
            "Uncertainties and errors"     => "man/uncertainties.md",
            "Properties calculations"      => "man/calculations.md",
            "Energy-resolved yields"       => "man/yields.md",
            "Chimeric spectra"             => "man/chimeric.md",
            "Interoperability"             => "man/interop.md",
            "Plotting"                     => "man/plotting.md",
        ],
        "References"                       => "reference.md",
        "Miscellaneous"                    => "misc.md",
    ],
)

deploydocs(
    repo = "github.com/ajgiuliani/MassJ.jl.git",
    devbranch = "master",
    push_preview = true,
)
