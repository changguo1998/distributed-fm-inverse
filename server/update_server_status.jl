#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)

if !isfile(FLAG_SERVER_UPLOADED())
    release_single_process_lock(@__FILE__)
    exit(0)
end

result_can_be_downloaded = filter(readdir(BUFFER_SERVER_RESULT())) do e
    tag = replace(e, "_result.tar.gz"=>"")
    edir = joinpath(BUFFER_SERVER_RUN(), tag)
    return isfile(joinpath(edir, FLAG_SERVER_PACK_RESULT_END))
end

info = Dict(
    "input" => map(t->replace(t, "_input.tar.gz"=>""), readdir(BUFFER_SERVER_INPUT())),
    "inverse" => readdir(BUFFER_SERVER_RUN()),
    "result" => map(t->replace(t, "_result.tar.gz"=>""), result_can_be_downloaded),
    "update_time" => now()
)

open(io->TOML.print(io, info), STATUS_SERVER(), "w")
release_single_process_lock(@__FILE__)