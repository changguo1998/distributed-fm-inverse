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

mkpath(BUFFER_HOST_UPLOAD)
mkpath(BUFFER_HOST_RESULT)
mkpath(BUFFER_SERVER_INPUT())
mkpath(BUFFER_SERVER_RUN())
mkpath(BUFFER_SERVER_RESULT())
mkpath(joinpath(PRJ_ROOT_PATH, "log"))
touch(FLAG_SERVER_UPLOADED())
