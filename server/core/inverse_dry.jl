using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using DelimitedFiles, SeisTools, TOML, Dates, JuliaSourceMechanism,
    Printf, JLD2, Statistics, LinearAlgebra

Geo = SeisTools.Geodesy

include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "multistage_lib.jl"))
include(joinpath(@__DIR__, "NLLOCshell.jl"))
include(joinpath(@__DIR__, "glibio.jl"))

edir = ARGS[1]

mkpath(joinpath(edir, "result"))

open(io->println(io, "Hello World!"), joinpath(edir, "result", "test.txt"), "w")

sleep(rand(5:10));

