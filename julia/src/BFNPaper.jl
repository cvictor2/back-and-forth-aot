module BFNPaper
using FFTW, MAT, JSON3, LinearAlgebra, Printf, SHA, Dates
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include("config.jl")
include("models.jl")
include("solver.jl")
include("io.jl")
include("plotting.jl")
include("suites.jl")
export Config, preset, run_experiment, run_suite, replot_paper, replot_result
end
