#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
const LOG_SETTING = (log=LOG_SERVER_INVERSE, lock=LOCK_SERVER_INVERSE_LOG)

buffer_run_list = readdir(BUFFER_SERVER_RUN)

if isempty(buffer_run_list)
    exit(0)
end

wait_for_submit = filter(buffer_run_list) do d
    edir = joinpath(BUFFER_SERVER_RUN, d)
    return isfile(joinpath(edir, FLAG_SERVER_UNPACKED)) && (!isfile(joinpath(edir, FLAG_SERVER_INVERSE_BEGIN)))
end

if isempty(wait_for_submit)
    exit(0)
end

evt = first(wait_for_submit)

log_info("call inverse for event $evt")

touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_BEGIN))

svrsetting = TOML.parsefile(SERVER_SETTING_FILE)

if DRY_RUN
    cmd = Cmd(`mkdir -p $(joinpath(BUFFER_SERVER_RUN, evt, "result")); sleep 5`)
else
    cmd = Cmd(["julia", "-t", svrsetting["threads_per_event"], joinpath(@__DIR__, "inverse.jl"), joinpath(BUFFER_SERVER_RUN, evt)])
end

try
    run(cmd)
    touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_END))
    log_info("inversion for event $evt succeeded")
catch err
    touch(joinpath(BUFFER_SERVER_RUN, evt, FLAG_SERVER_INVERSE_FAILED))
    log_error("inversion for event $evt failed")
    error(err)
end
