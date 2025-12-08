#!/usr/bin/env julia

# usage
# julia upload_to_server.jl

include(joinpath(@__DIR__, "../lib.jl"))
get_single_process_lock(@__FILE__)
const LOG_SETTING = (log=LOG_HOST_UPLOAD, lock=LOCK_HOST_UPLOAD_LOG)

buffered_event = readdir(BUFFER_HOST_UPLOAD)

if isempty(buffered_event)
    release_single_process_lock(@__FILE__)
    exit(0)
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

if maximum(priority) < 0
    release_single_process_lock(@__FILE__)
    exit(0)
end

i = argmax(priority)
svr = nodes.servers[i]
datafile = first(buffered_event)
tag = replace(datafile, "_input.tar.gz"=>"")
server_input_buffer = replace(BUFFER_SERVER_INPUT, PRJ_ROOT_PATH=>svr.system_root)
scpfrom = joinpath(BUFFER_HOST_UPLOAD, datafile)
scpto = joinpath(server_input_buffer, datafile)

if DEBUG
    @info "upload file from $scpfrom to $scpto"
end

if svr.hostname == nodes.host.hostname
    cmd1 = Cmd(["rm", FLAG_SERVER_UPLOADED])
    cmd2 = Cmd(["cp", scpfrom, scpto])
    cmd3 = Cmd(["touch",  FLAG_SERVER_UPLOADED])
else
    server_address = svr.user*"@"*svr.ip
    server_upload_flag_file = replace(FLAG_SERVER_UPLOADED, PRJ_ROOT_PATH=>svr.system_root)

    cmd1 = Cmd(["ssh", server_address, "rm \"$(server_upload_flag_file)\""])
    cmd2 = Cmd(["scp", scpfrom, "$server_address:$scpto"])
    cmd3 = Cmd(["ssh", server_address, "touch \"$(server_upload_flag_file)\""])
end

try
    run(cmd1)
    log_info("update upload flag: ", string(cmd1))
    run(cmd2)
    log_info("send data: ", string(cmd2))
    rm(scpfrom; force=true)
    run(cmd3)
    log_info("update upload flag: ", string(cmd3))
    log_info("upload data file $datafile to server $(svr.hostname)")
catch err
    log_err("failed to send data to server $(svr.hostname)")
    error(err)
finally
    release_single_process_lock(@__FILE__)
end

