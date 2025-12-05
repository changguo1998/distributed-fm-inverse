#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
randstr(n::Integer=4) = join(rand(['A':'Z'; 'a':'z'; '0':'9'], n))
N = isempty(ARGS) ? 10 : parse(Int, ARGS[1])

for i = 1:N
    edir = joinpath(@__DIR__, "monitored_dir", "TestEvent$(i)_$(randstr())")
    mkpath(joinpath(edir, "greenfun"))
    touch(joinpath(edir, FLAG_HOST_PREPROCESS_BEGIN))
    touch(joinpath(edir, "auto.jld2"))
    touch(joinpath(edir, FLAG_HOST_PREPROCESS_END))
end
