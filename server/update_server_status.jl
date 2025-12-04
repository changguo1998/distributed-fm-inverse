#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

get_single_process_lock(@__DIR__)

while true
    if isfile(FLAG_SERVER_UPLOADED)
        break
    end
    sleep(rand(1:5))
end

open(STATUS_SERVER, "w") do io
    TOML.print(io, Dict(
        "input" => map(t->replace(t, "_input.tar.gz"=>""), readdir(BUFFER_SERVER_INPUT)),
        "inverse" => readdir(BUFFER_SERVER_RUN),
        "result" => map(t->replace(t, "_result.tar.gz"=>""), readdir(BUFFER_SERVER_RESULT)),
        "update_time" => now()
    ))
end

release_single_process_lock(@__DIR__)