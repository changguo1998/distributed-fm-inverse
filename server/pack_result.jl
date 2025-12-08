#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)
const LOG_SETTING = (log=LOG_SERVER_RESULT, lock=LOCK_SERVER_RESULT_LOG)

event_finished = filter(readdir(BUFFER_SERVER_RUN())) do e
    return isfile(joinpath(BUFFER_SERVER_RUN(), e, FLAG_SERVER_INVERSE_END)) ||
        isfile(joinpath(BUFFER_SERVER_RUN(), e, FLAG_SERVER_INVERSE_FAILED))
end

if isempty(event_finished)
    release_single_process_lock(@__FILE__)
    exit(0)
end

event_already_packed = map(f->replace(f, "_result.tar.gz"=>""), readdir(BUFFER_SERVER_RESULT()))
wait_for_pack = setdiff(event_finished, event_already_packed)

if isempty(wait_for_pack)
    release_single_process_lock(@__FILE__)
    exit(0)
end

evt = first(wait_for_pack)
edir = joinpath(BUFFER_SERVER_RUN(), evt)

touch(joinpath(edir, FLAG_SERVER_PACK_RESULT_BEGIN))

if !isdir(joinpath(edir, "result"))
    touch(joinpath(edir, FLAG_SERVER_INVERSE_FAILED))
end

if isempty(readdir(joinpath(edir, "result")))
    touch(joinpath(edir, FLAG_SERVER_INVERSE_FAILED))
end

if isfile(joinpath(edir, FLAG_SERVER_INVERSE_FAILED))
    touch(joinpath(edir, "inversion_failed.flag"))
    cmd = Cmd(Cmd(["tar", "czf", joinpath(BUFFER_SERVER_RESULT(), evt*"_result.tar.gz"), "inversion_failed.flag"]); dir=edir)
else
    cmd = Cmd(Cmd(["tar", "czf", joinpath(BUFFER_SERVER_RESULT(), evt*"_result.tar.gz"), "result"]); dir=edir)
end

try
    run(cmd)
    touch(joinpath(edir, FLAG_SERVER_PACK_RESULT_END))
    log_info("pack result for event $evt")
catch err
    log_err("failed to pack result for event $evt")
    error(err)
finally
    release_single_process_lock(@__FILE__)
end
