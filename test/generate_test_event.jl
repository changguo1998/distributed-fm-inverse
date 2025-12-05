#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

N = 10

for i = 1:N
    mkpath(joinpath(@__DIR__, "monitored_dir", "TestEvent$i"))
    touch(joinpath(@__DIR__, "monitored_dir", "TestEvent$i", "test_file.txt"))
end