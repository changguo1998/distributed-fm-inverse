#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

const LOG_SETTING = (log=LOG_SERVER_RESULT, lock=LOCK_SERVER_RESULT_LOG)

event_finished = readdir(BUFFER_SERVER_RUN)

if isempty(event_finished)
    exit(0)
end

event_already_packed = map(f->replace(f, "_result.tar.gz"=>""), readdir(BUFFER_SERVER_RESULT))
wait_for_pack = setdiff(event_finished, event_already_packed)

if isempty(wait_for_pack)
    exit(0)
end

evt = first(wait_for_pack)

log_info("pack result for event $evt")

edir = joinpath(BUFFER_SERVER_RUN, evt)

touch(joinpath(edir, FLAG_SERVER_RESULT_BEGIN))

if isfile(joinpath(edir, FLAG_SERVER_INVERSE_FAILED))
    touch(joinpath(edir, "inversion_failed.flag"))
    cmd = Cmd(Cmd(["tar", "czf", joinpath(BUFFER_SERVER_RESULT, evt*"_result.tar.gz"), "inversion_failed.flag"]); dir=edir)
else
    cmd = Cmd(Cmd(["tar", "czf", joinpath(BUFFER_SERVER_RESULT, evt*"_result.tar.gz"), "result"]); dir=edir)
end

try
    run(cmd)
catch
    log_error("failed to pack result for event $evt")
    exit(0)
end

touch(joinpath(edir, FLAG_SERVER_RESULT_END))
