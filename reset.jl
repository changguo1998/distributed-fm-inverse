#!/usr/bin/env julia
include(joinpath(@__DIR__, "lib.jl"))

function cleanup(dir::AbstractString)
    if isdir(joinpath(PRJ_ROOT_PATH, dir))
        @info dir
        foreach(f -> rm(f; recursive=true),
            readdir(joinpath(PRJ_ROOT_PATH, dir); join=true))
    end
end

cleanup("log")
cleanup("var")

run(`julia $(joinpath(PRJ_ROOT_PATH, "install.jl"))`, devnull, devnull, devnull)

touch(FLAG_SERVER_UPLOADED)