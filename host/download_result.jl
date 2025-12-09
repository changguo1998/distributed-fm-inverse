#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
LOG_SETTING = (log=LOG_HOST_DOWNLOAD, lock=LOCK_HOST_DOWNLOAD_LOG)

function main()
    nodes = host_load_node()

    for svr in nodes.servers
        status = get_server_loading(svr, nodes.host)
        if isnothing(status)
            continue
        end
        if isempty(status["result"])
            continue
        end
        svr_address = svr.user * "@" * svr.ip
        svr_input_buffer = BUFFER_SERVER_INPUT(svr)
        svr_run_buffer = BUFFER_SERVER_RUN(svr)
        svr_result_buffer = BUFFER_SERVER_RESULT(svr)
        svr_clean_buffer = BUFFER_SERVER_CLEAN(svr)
        for r in status["result"]
            if svr.hostname == nodes.host.hostname
                cmd_scp = Cmd(["cp", joinpath(BUFFER_SERVER_RESULT(), r*"_result.tar.gz"), BUFFER_HOST_RESULT])
                cmd_rm = Cmd(["touch", joinpath(svr_clean_buffer, r * ".flag")])
            else
                cmd_scp = Cmd(["scp", svr_address*":"*abspath(svr_result_buffer, r*"_result.tar.gz"), BUFFER_HOST_RESULT])
                cmd_rm = Cmd(["ssh", svr_address, "touch $(joinpath(svr_clean_buffer, r * ".flag"))"])
            end
            if DEBUG
                @info cmd_scp
                @info cmd_rm
            end
            try
                run(cmd_scp)
                log_info("download from $(svr.hostname)")
                run(cmd_rm)
            catch err
                log_err("error while downloading from $(svr.hostname)")
            end
        end
    end

    host_wait_for_unpack = readdir(BUFFER_HOST_RESULT)

    if isempty(host_wait_for_unpack)
        return nothing
    end

    queue_info = nothing
    try
        get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
        queue_info = TOML.parsefile(STATUS_QUEUE)
    catch err
        error(err)
    finally
        release_lock(LOCK_HOST_QUEUE_STATUS_FILE)
    end
    key_remove = String[]
    for f in host_wait_for_unpack
        tag = replace(f, "_result.tar.gz"=>"")
        if !haskey(queue_info, tag)
            log_err("Cannot find $tag in queue info")
            rm(joinpath(BUFFER_HOST_RESULT, f); force=true)
            continue
        end
        cmd = Cmd([
            "tar",
            "xaf",
            joinpath(BUFFER_HOST_RESULT, f),
            "-C",
            queue_info[tag]
        ])
        try
            run(cmd)
            touch(joinpath(queue_info[tag], FLAG_HOST_INVERSION_FINISHED))
            push!(key_remove, tag)
            rm(joinpath(BUFFER_HOST_RESULT, f); force=true)
            log_info("unpack $tag to $(queue_info[tag])")
        catch
            log_err("error while unpacking $f")
        end
    end

    # get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
    # q = TOML.parsefile(STATUS_QUEUE)
    # t = Dict()
    # for k in keys(q)
    #     if k in key_remove
    #         continue
    #     end
    #     t[k] = q[k]
    # end
    # t["update_time"] = now()
    # open(io->TOML.print(io, t), STATUS_QUEUE, "w")
    # release_lock(LOCK_HOST_QUEUE_STATUS_FILE)
end

get_single_process_lock(@__FILE__)
try
    main()
catch err
    log_err("failed to run script")
    log_err(string(err))
end
release_single_process_lock(@__FILE__)