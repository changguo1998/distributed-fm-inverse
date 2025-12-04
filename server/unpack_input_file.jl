#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

const LOG_SETTING = (log=LOG_SERVER_UNPACK, lock=LOCK_SERVER_UNPACK_LOG)

if !isfile(FLAG_SERVER_UPLOADED)
    exit(0)
end

buffer_list = readdir(BUFFER_SERVER_INPUT)
inv_list = readdir(BUFFER_SERVER_RUN)
wait_for_submit = setdiff(buffer_list, inv_list)

if isempty(wait_for_submit)
    exit(0)
end

for f in wait_for_submit
    log_info("unpack $f")
    tag = replace(f, "_input.tar.gz"=>"")
    run_dir = joinpath(BUFFER_SERVER_RUN, tag)
    mkpath(run_dir)
    cmd = Cmd(["tar", "xaf", joinpath(BUFFER_SERVER_INPUT, f), "-C", run_dir])
    try
        run(cmd)
        touch(joinpath(run_dir, FLAG_SERVER_UNPACKED))
        log_info("$f unpacked successfully")
    catch
        log_error("failed to unpack $f")
        continue
    end
end