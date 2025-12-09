#!/usr/bin/env julia
include(joinpath(@__DIR__, "lib.jl"))

using ArgumentProcessor, CairoMakie, DWN, Dates, DelimitedFiles, JLD2,
    JuliaSourceMechanism, LengthAreaVolume, LinearAlgebra, Printf, SHA,
    SeisTools, SeismicRayTrace, Statistics, TOML

nodes = host_load_node()

