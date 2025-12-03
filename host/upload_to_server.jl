#!/usr/bin/env julia

# usage
# julia upload_to_server.jl

include(joinpath(@__DIR__, "../lib.jl"))

const LOG_SETTING = (log=LOG_HOST_UPLOAD, lock=LOCK_HOST_UPLOAD_LOG)

log_info("upload to server start")

buffered_event = readdir(BUFFER_HOST_UPLOAD)

if isempty(buffered_event)
    log_info("no event to upload, exit")
    exit(0)
end

nodes = host_load_node()
priority = map(nodes.servers) do svr
    status = get_server_loading(svr)
    if isnothing(status)
        return -1
    end
    remaining = svr.max_event_number - length(status.input)
    if remaining < 0.5
        return -1
    end
    return svr.priority * 10 + remaining
end

if maximum(priority) < 0
    log_info("server loading is full, exit")
    exit(0)
end

i = argmax(priority)
svr = nodes.servers[i]
datafile = first(buffered_event)
tag = replace(datafile, "_input.tar.gz"=>"")

server_address = svr.user*"@"*svr.ip
server_upload_flag_file = replace(FLAG_SERVER_UPLOADED, PRJ_ROOT_PATH=>svr.system_root)
server_input_buffer = replace(BUFFER_SERVER_INPUT, PRJ_ROOT_PATH=>svr.system_root)
scpfrom = joinpath(BUFFER_HOST_UPLOAD, datafile)
scpto = joinpath(server_input_buffer, datafile)

log_info("upload data file $datafile to server $(svr.hostname)")
try
    run(`ssh $(server_address) "rm $(server_upload_flag_file)"`)
    run(`scp $scpfrom "$server_address:$scpto"`)
    rm(scpfrom; force=true)
    run(`ssh $(server_address) "touch $(server_upload_flag_file)"`)
catch err
    log_err("failed to send data to server $(svr.hostname)")
    error(err)
end

log_info("upload to server done")
