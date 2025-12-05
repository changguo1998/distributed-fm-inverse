#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib.jl"))

svr_setting = TOML.parsefile(SERVER_SETTING_FILE)
nodes = host_load_node()

planed_script = Tuple{String,Second}[]

if svr_setting["hostname"] == nodes.host.hostname
    append!(planed_script, [
        ("host/submit_to_host_buffer.jl", Second(30)),
        ("host/upload_to_server.jl", Second(30)),
        ("host/download_result.jl", Second(30)),
    ])
end

append!(planed_script, [
    ("server/unpack_input_file.jl", Second(30)),
    ("server/call_inverse.jl", Second(30)),
    ("server/pack_result.jl", Second(30)),
    ("server/update_server_status.jl", Second(10))
])

lastrun = fill(now() - Second(1) - maximum(getindex.(planed_script, 2)), length(planed_script))

while true
    global lastrun
    local t = now()
    for i = eachindex(planed_script)
        if t - lastrun[i] >= planed_script[i][2]
            @info "[$t] run script: $(planed_script[i][1])"
            cmd = Cmd(`julia $(planed_script[i][1])`; dir=@__DIR__)
            run(cmd, devnull, devnull, devnull; wait=false)
            lastrun[i] = t
        end
    end
    sleep(1)
end