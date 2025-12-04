#!/usr/bin/env julia

# usage
# julia submit_to_host_buffer.jl /absolute/path/to/event

include(joinpath(@__DIR__, "../lib.jl"))

get_single_process_lock(@__DIR__)

const LOG_SETTING = (log=LOG_HOST_QUEUE, lock=LOCK_HOST_QUEUE_LOG)

randstr(n::Integer) = join(rand(['A':'Z', 'a':'z', '0':'9'], n), "")

edir = abspath(ARGS[1])

# check input data status
if !isfile(joinpath(edir, FLAG_HOST_PREPROCESS_END))
    log_warn("event $edir not processed")
    release_single_process_lock(@__DIR__)
    exit(0)
end
if !isfile(joinpath(edir, "auto.jld2"))
    log_err("event $edir cannot find auto.jld2")
    release_single_process_lock(@__DIR__)
    error("data not ready")
end
if !isdir(joinpath(edir, "greenfun"))
    log_err("event $edir cannot find greenfun")
    release_single_process_lock(@__DIR__)
    error("data not ready")
end

if isfile(joinpath(edir, FLAG_HOST_QUEUE_BEGIN))
    log_warn("event $edir queue already begins")
    release_single_process_lock(@__DIR__)
    exit(0)
end

touch(joinpath(edir, FLAG_HOST_QUEUE_BEGIN))

tag = ""
while true
    global tag
    tag = randstr(16)
    if !isfile(abspath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz"))
        break
    end
end

log_info("tag: $tag")

cmd = Cmd(`tar czf $(abspath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz")) auto.jld2 greenfun`; dir=edir)

try
    run(cmd)
    log_info("pack up ", tag, " ", string(cmd))
catch err
    log_err("pack up failed")
    release_single_process_lock(@__DIR__)
    error(err)
end

get_lock(LOCK_HOST_QUEUE_STATUS_FILE)

if isfile(STATUS_QUEUE)
    t = TOML.parsefile(STATUS_QUEUE)
else
    t = Dict()
end

t[tag] = edir
t["update_time"] = now()

open(io->TOML.print(io, t), STATUS_QUEUE, "w")
release_lock(LOCK_HOST_QUEUE_STATUS_FILE)

touch(joinpath(edir, FLAG_HOST_QUEUE_END))

log_info("submit $edir to queue")
release_single_process_lock(@__DIR__)
