#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
const LOG_SETTING = (log=LOG_HOST_UPLOAD, lock=LOCK_HOST_UPLOAD_LOG)
function main()
    if !isfile(STATUS_QUEUE)
        if DEBUG
            @info "No status queue file found. Exiting."
        end
        return nothing
    end

    t = nothing
    try
        t = TOML.parsefile(STATUS_QUEUE)
    catch err
        log_err("error while reading STATUS_QUEUE file.")
        error(err)
    end

    if isnothing(t)
        if DEBUG
            @info "Failed to parse STATUS_QUEUE file. Exiting."
        end
        return nothing
    end

    if length(collect(keys(t))) < 2
        if DEBUG
            @info "No events in queue status. Exiting."
        end
        return nothing
    end
    event_zipped = collect(keys(t))
    event_not_uploaded = map(f->replace(f, "_input.tar.gz"=>""), readdir(BUFFER_HOST_UPLOAD))
    buffered_event = intersect(event_zipped, event_not_uploaded)

    if isempty(buffered_event)
        if DEBUG
            @info "No events to upload. Exiting."
        end
        return nothing
    end

    tag = first(buffered_event)

    if DEBUG
        @info "Uploading $tag"
        @info "Check server status"
    end

    nodes = host_load_node()
    priority = map(nodes.servers) do svr
        status = get_server_loading(svr, nodes.host)
        if isnothing(status)
            return -1
        end
        remaining = svr.max_event_number - length(status["input"])
        if remaining < 0.5
            return -1
        end
        return svr.priority * 100 + remaining
    end

    if DEBUG
        @info "Priority: $priority"
    end

    if maximum(priority) < 0
        if DEBUG
            @info "No available server. Exiting."
        end
        return nothing
    end

    i = argmax(priority)
    svr = nodes.servers[i]

    if DEBUG
        @info "Selected server: $svr"
    end

    scpfrom = joinpath(BUFFER_HOST_UPLOAD, tag*"_input.tar.gz")
    scpto = joinpath(BUFFER_SERVER_INPUT(svr), tag*"_input.tar.gz")

    if DEBUG
        @info "upload file from $scpfrom to $scpto"
    end

    if svr.hostname == nodes.host.hostname
        cmd1 = Cmd(["rm", FLAG_SERVER_UPLOADED()])
        cmd2 = Cmd(["cp", scpfrom, scpto])
        cmd3 = Cmd(["touch",  FLAG_SERVER_UPLOADED()])
    else
        server_address = svr.user*"@"*svr.ip
        server_upload_flag_file = FLAG_SERVER_UPLOADED(svr)

        cmd1 = Cmd(["ssh", server_address, "rm \"$(server_upload_flag_file)\""])
        cmd2 = Cmd(["scp", scpfrom, "$server_address:$scpto"])
        cmd3 = Cmd(["ssh", server_address, "touch \"$(server_upload_flag_file)\""])
    end

    try
        run(cmd1)
        if DEBUG
            log_info("update upload flag: ", string(cmd1))
        end
    catch
        log_warn("rm server $(svr.hostname) upload flag failed")
    end

    try
        run(cmd2, devnull, devnull, devnull)
        if DEBUG
            log_info("send data: ", string(cmd2))
        end
        rm(scpfrom; force=true)
        log_info("upload data file $tag to server $(svr.hostname)")
    catch err
        log_err("failed to send data to server $(svr.hostname)")
        error(err)
    end

    try
        run(cmd3)
        if DEBUG
            log_info("update upload flag: ", string(cmd3))
        end
    catch
        log_warn("touch server $(svr.hostname) upload flag failed")
    end
end


try
    get_single_process_lock(@__FILE__)
    main()
catch err
    log_err("failed to run script")
    log_err(string(err))
    error(err)
finally
    release_single_process_lock(@__FILE__)
end