#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

LOG_SETTING = (log=LOG_HOST_DOWNLOAD, lock=LOCK_HOST_DOWNLOAD_LOG)

nodes = host_load_node()

rm_cmds = Cmd[]

for svr in nodes.servers
    status = get_server_loading(svr)
    if isnothing(status)
        continue
    end
    if isempty(status["result"])
        continue
    end
    svr_address = svr.user * "@" * svr.ip
    svr_result_buffer = replace(BUFFER_SERVER_RESULT, PRJ_ROOT_PATH=>svr.system_root)
    cmd_scp = Cmd([
        "scp";
        map(r->svr_address*":"*joinpath(svr_result_buffer, r*"_result.tar.gz"), status["result"]);
        BUFFER_HOST_RESULT
    ])
    cmd_rm = Cmd([
        "ssh",
        svr_address,
        """
rm $(map(r->joinpath(svr_result_buffer, r*"_result.tar.gz"), status["result"]))
        """
    ])
    push!(rm_cmds, cmd_rm)
    log_info("download from $(svr.hostname)")
    try
        run(cmd)
    catch err
        log_err("error while downloading from $(svr.hostname)")
    end
end

log_info("unpack result")
get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
queue_info = TOML.parsefile(STATUS_QUEUE)
release_lock(LOCK_HOST_QUEUE_STATUS_FILE)
key_remove = String[]
for f in readdir(BUFFER_HOST_RESULT)
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
    log_info("unpack $tag to $(queue_info[tag])")
    try
        run(cmd)
    catch
        log_err("error while unpacking $f")
    end
    push!(key_remove, tag)
end

log_info("update queue status")
get_lock(LOCK_HOST_QUEUE_STATUS_FILE)
q = TOML.parsefile(STATUS_QUEUE)
t = Dict()
for k in keys(q)
    if k in key_remove
        continue
    end
    t[k] = q[k]
end
open(io->TOML.print(io, t), STATUS_QUEUE, "w")
release_lock(LOCK_HOST_QUEUE_STATUS_FILE)

log_info("download result done")