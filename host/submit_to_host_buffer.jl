#!/usr/bin/env julia

# usage
# julia submit_to_host_buffer.jl /absolute/path/to/event

include(joinpath(@__DIR__, "../lib.jl"))
randstr(n::Integer) = join(rand(['A':'Z'; 'a':'z'; '0':'9'], n))
get_single_process_lock(@__FILE__)
const LOG_SETTING = (log=LOG_HOST_QUEUE, lock=LOCK_HOST_QUEUE_LOG)

nodes = host_load_node()
candidate_events = String[]
for d in nodes.host.monitor_directory
    if !isdir(d)
        continue
    end
    tmp_list = filter(readdir(d; join=true)) do e
        return isfile(joinpath(e, FLAG_HOST_PREPROCESS_END)) &&
           isfile(joinpath(e, "auto.jld2")) &&
           isdir(joinpath(e, "greenfun")) &&
           (!isfile(joinpath(e, FLAG_HOST_QUEUE_BEGIN))) &&
           (!isfile(joinpath(e, FLAG_HOST_INVERSION_FINISHED)))
    end
    if isempty(tmp_list)
        continue
    end
    append!(candidate_events, tmp_list)
end

if isempty(candidate_events)
    release_single_process_lock(@__FILE__)
    exit(0)
end

edir = rand(candidate_events)

touch(joinpath(edir, FLAG_HOST_QUEUE_BEGIN))

tag = ""
while true
    global tag
    tag = randstr(16)
    if !isfile(abspath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz"))
        break
    end
end

cmd = Cmd(`tar czf $(abspath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz")) auto.jld2 greenfun`; dir=edir)

try
    run(cmd)
    log_info("enqueue ", tag, ", from ", edir)
catch err
    log_err("failed to enqueue ", tag, ", from ", edir)
    release_single_process_lock(@__FILE__)
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
release_single_process_lock(@__FILE__)
