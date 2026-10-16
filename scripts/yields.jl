#!/usr/bin/env julia
# yields.jl — CLI wrapper around MassJ.yields for energy-resolved peak integration.
#
# For every mass-spectrum file found in --input, the script:
#   1. averages the scans into one spectrum,
#   2. integrates each (mz1, mz2) window defined in the peak list,
#   3. sums all peak integrals into a row-TIC,
#   4. assigns an external parameter value x = x0 + step * i (i = 0..N-1,
#      natural-sorted files),
#   5. writes the table as CSV.
#
# Usage:
#     julia --project=@. scripts/yields.jl --input <dir> --list peaks.csv \
#           --energy E0 --step st [--output out.csv] [--tic] [--flux flux.txt]
#
# The peak list CSV is read by MassJ.read_peaklist: 2 columns (mz, label), 3 columns
# (mz1, mz2, label) or 4 columns (mz, tol, method, label), with or without header.

using MassJ

const HELP = """
yields.jl — CLI for MassJ.yields (energy-resolved peak integration).

Required arguments:
  --input   <dir>    directory containing one spectrum file per energy step
  --list    <csv>    peak list CSV: (mz, label), (mz1, mz2, label) or
                     (mz, tol, method, label)
  --energy  <float>  initial value of the external parameter (x0)
  --step    <float>  increment between successive files

Optional:
  --output  <csv>    output CSV path (default: <input>/integrated_yields.csv)
  --xlabel  <str>    x-axis label written to the CSV (default: "energy")
  --tic              normalize peak columns by row-TIC (peaks sum to 1; raw TIC kept)
  --flux    <txt>    divide every column by the photon flux, linearly interpolated at
                     the row's x value (MassJ.normalize_external). Flux file: header
                     lines and # comments skipped; columns x, flux[, 1-sigma].
  --help             show this message
"""

function parse_args(args::Vector{String})
    opts = Dict{String,Any}("tic" => false, "xlabel" => "energy")
    i = 1
    while i <= length(args)
        a = args[i]
        if     a in ("--input",  "-i"); opts["input"]  = args[i+1]; i += 2
        elseif a in ("--list",   "-l"); opts["list"]   = args[i+1]; i += 2
        elseif a in ("--energy", "-e"); opts["energy"] = parse(Float64, args[i+1]); i += 2
        elseif a in ("--step",   "-s"); opts["step"]   = parse(Float64, args[i+1]); i += 2
        elseif a in ("--output", "-o"); opts["output"] = args[i+1]; i += 2
        elseif a in ("--flux",   "-f"); opts["flux"]   = args[i+1]; i += 2
        elseif a == "--xlabel";         opts["xlabel"] = args[i+1]; i += 2
        elseif a in ("--tic",    "-t"); opts["tic"]    = true;      i += 1
        elseif a in ("--help",   "-h"); print(HELP); exit(0)
        else
            println(stderr, "Unknown argument: $a")
            print(stderr, HELP); exit(2)
        end
    end
    for k in ("input", "list", "energy", "step")
        haskey(opts, k) || (println(stderr, "Missing --$k"); print(stderr, HELP); exit(2))
    end
    return opts
end

function main()
    opts = parse_args(ARGS)

    peaks = MassJ.read_peaklist(opts["list"])
    yc = MassJ.yields(opts["input"], peaks;
                      x0     = opts["energy"],
                      step   = opts["step"],
                      xlabel = opts["xlabel"])
    @info "Integrated $(length(yc.x)) files × $(length(yc.labels)) peaks"

    if opts["tic"]
        yc = MassJ.normalize_tic(yc)
        @info "Normalized peak columns by row-TIC"
    end

    if haskey(opts, "flux")
        yc = MassJ.normalize_external(yc, opts["flux"])
        @info "Flux-normalized using $(opts["flux"])"
    end

    outpath = get(opts, "output", joinpath(opts["input"], "integrated_yields.csv"))
    MassJ.write_csv(yc, outpath)
    @info "Wrote $outpath  ($(length(yc.x)) rows × $(2 * length(yc.labels) + 3) cols, with 1-σ columns)"
end

main()
