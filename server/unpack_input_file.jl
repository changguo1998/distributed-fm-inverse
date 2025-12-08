#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)
const LOG_SETTING = (log=LOG_SERVER_UNPACK, lock=LOCK_SERVER_UNPACK_LOG)

if !isfile(FLAG_SERVER_UPLOADED())
    release_single_process_lock(@__FILE__)
    exit(0)
end

buffer_list = map(t->replace(t, "_input.tar.gz"=>""), readdir(BUFFER_SERVER_INPUT()))
inv_list = readdir(BUFFER_SERVER_RUN())
wait_for_submit = setdiff(buffer_list, inv_list)

if isempty(wait_for_submit)
    release_single_process_lock(@__FILE__)
    exit(0)
end

tag = first(wait_for_submit)

log_info("unpack $tag")
run_dir = joinpath(BUFFER_SERVER_RUN(), tag)
mkpath(run_dir)
cmd = Cmd(["tar", "xaf", joinpath(BUFFER_SERVER_INPUT(), tag * "_input.tar.gz"), "-C", run_dir])
try
    run(cmd)
    touch(joinpath(run_dir, FLAG_SERVER_UNPACKED))
    log_info("$tag unpacked successfully")
catch e
    log_err("failed to unpack $tag")
    error(e)
finally
    release_single_process_lock(@__FILE__)
end