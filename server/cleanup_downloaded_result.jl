#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)

function _rm(f::AbstractString)
    rm(f; recursive=true, force=true)
    return nothing
end

rm_list = filter(readdir(BUFFER_SERVER_CLEAN())) do f
    s = stat(joinpath(BUFFER_SERVER_CLEAN(), f))
    return now() - unix2datetime(max(s.mtime, s.ctime)) > Second(2)
end

foreach(rm_list) do f
    tag = replace(f, ".flag"=>"")
    _rm(joinpath(BUFFER_SERVER_INPUT(), tag * "_input.tar.gz"))
    _rm(joinpath(BUFFER_SERVER_RUN(), tag))
    _rm(joinpath(BUFFER_SERVER_RESULT(), tag * "_result.tar.gz"))
    _rm(joinpath(BUFFER_SERVER_CLEAN(), f))
end

release_single_process_lock(@__FILE__)