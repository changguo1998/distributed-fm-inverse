#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
randstr(n::Integer) = join(rand(['A':'Z'; 'a':'z'; '0':'9'], n))
const LOG_SETTING = (log=LOG_HOST_QUEUE, lock=LOCK_HOST_QUEUE_LOG)
function main()
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
            isdir(joinpath(e, "sac")) &&
            (!isfile(joinpath(e, FLAG_HOST_QUEUE_BEGIN))) &&
            (!isfile(joinpath(e, FLAG_HOST_INVERSION_FINISHED))) &&
            (!isdir(joinpath(e, "result")))
        end
        if isempty(tmp_list)
            continue
        end
        append!(candidate_events, tmp_list)
    end

    if isempty(candidate_events)
        return nothing
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

    cmd = Cmd(`tar czf $(abspath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz")) auto.jld2 greenfun sac`; dir=edir)

    try
        run(cmd)
        log_info("enqueue ", tag, ", from ", edir)
    catch err
        log_err("failed to enqueue ", tag, ", from ", edir)
        error(err)
    end

    t = Dict()

    if isfile(STATUS_QUEUE)
        try
            get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
            t = TOML.parsefile(STATUS_QUEUE)
        catch err
            release_lock(LOCK_HOST_QUEUE_STATUS_FILE)
            log_err("failed to read ", STATUS_QUEUE)
            error(err)
        end
    end

    t[tag] = edir
    t["update_time"] = now()

    get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
    open(io->TOML.print(io, t), STATUS_QUEUE, "w")
    release_lock(LOCK_HOST_QUEUE_STATUS_FILE)

    touch(joinpath(edir, FLAG_HOST_QUEUE_END))
end

get_single_process_lock(@__FILE__)
get_lock(LOCK_HOST_STATUS_UPLOADING)
try
    main()
catch err
    log_err("failed to run script")
    log_err(string(err))
end
release_lock(LOCK_HOST_STATUS_UPLOADING)
release_single_process_lock(@__FILE__)
