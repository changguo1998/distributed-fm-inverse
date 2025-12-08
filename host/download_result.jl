#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)
LOG_SETTING = (log=LOG_HOST_DOWNLOAD, lock=LOCK_HOST_DOWNLOAD_LOG)

nodes = host_load_node()

rm_cmds = Cmd[]

for svr in nodes.servers
    status = get_server_loading(svr, nodes.host)
    if isnothing(status)
        continue
    end
    if isempty(status["result"])
        continue
    end
    if svr.hostname == nodes.host.hostname
        cmd_scp = Cmd([
            "cp",
            map(r->joinpath(BUFFER_SERVER_RESULT, r*"_result.tar.gz"), status["result"])...,
            BUFFER_HOST_RESULT
            ])
        cmd_rm = Cmd([
            "rm",
            map(r->joinpath(BUFFER_SERVER_RESULT, r*"_result.tar.gz"), status["result"])...
        ])
    else
        svr_address = svr.user * "@" * svr.ip
        svr_input_buffer = BUFFER_SERVER_INPUT(svr)
        svr_run_buffer = BUFFER_SERVER_RUN(svr)
        svr_result_buffer = BUFFER_SERVER_RESULT(svr)
        cmd_scp = Cmd([
            "scp";
            map(r->svr_address*":"*joinpath(svr_result_buffer, r*"_result.tar.gz"), status["result"])...;
            BUFFER_HOST_RESULT
        ])
        cmd_rm = Cmd([
            "ssh",
            svr_address,
            "rm $(map(r->joinpath(svr_input_buffer, r*"_input.tar.gz"), status["result"])...);",
            "rm -r $(status["result"]...);",
            "rm $(map(r->joinpath(svr_result_buffer, r*"_result.tar.gz"), status["result"])...);"
        ])
    end
    if DEBUG
        @info cmd_scp
        @info cmd_rm
    end
    try
        run(cmd_scp)
        push!(rm_cmds, cmd_rm)
        log_info("download from $(svr.hostname)")
    catch err
        log_err("error while downloading from $(svr.hostname)")
    end
end

foreach(rm_cmds) do c
    try
        run(c)
    finally
        return nothing
    end
end

host_wait_for_unpack = readdir(BUFFER_HOST_RESULT)

if isempty(host_wait_for_unpack)
    release_single_process_lock(@__FILE__)
    exit(0)
end

get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
queue_info = TOML.parsefile(STATUS_QUEUE)
release_lock(LOCK_HOST_QUEUE_STATUS_FILE)
key_remove = String[]
for f in host_wait_for_unpack
    tag = replace(f, "_result.tar.gz"=>"")
    if !haskey(queue_info, tag)
        log_err("Cannot find $tag in queue info")
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
        log_info("unpack $tag to $(queue_info[tag])")
    catch
        log_err("error while unpacking $f")
    end
end

get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
q = TOML.parsefile(STATUS_QUEUE)
t = Dict()
for k in keys(q)
    if k in key_remove
        continue
    end
    t[k] = q[k]
end
t["update_time"] = now()
open(io->TOML.print(io, t), STATUS_QUEUE, "w")
release_lock(LOCK_HOST_QUEUE_STATUS_FILE)

release_single_process_lock(@__FILE__)