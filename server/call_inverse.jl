#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

const LOG_SETTING = (log=LOG_SERVER_INVERSE, lock=LOCK_SERVER_INVERSE_LOG)

buffer_run_list = readdir(BUFFER_SERVER_RUN)

if isempty(buffer_run_list)
    exit(0)
end

wait_for_submit = filter(d->!isfile(joinpath(BUFFER_SERVER_RUN, d, FLAG_SERVER_UNPACKED)), buffer_run_list)

if isempty(wait_for_submit)
    exit(0)
end

evt = first(wait_for_submit)

log_info("call inverse for event $evt")

touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_BEGIN))

svrsetting = TOML.parsefile(SERVER_SETTING_FILE)

cmd = Cmd(["julia", "-t", svrsetting["threads_per_event"], joinpath(@__DIR__, "inverse.jl"), joinpath(BUFFER_SERVER_RUN, evt)])

try
    run(cmd)
catch
    log_error("inversion for event $evt failed")
    touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_FAILED))
    exit(0)
end

touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_END))