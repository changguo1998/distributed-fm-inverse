#!/usr/bin/env julia
include(joinpath(@__DIR__, "lib.jl"))

# reset
cmd = Cmd(`bash dfmi.sh reset.jl`; dir=@__DIR__)
if DEBUG
    @info cmd
end
run(cmd, devnull, devnull, devnull)

# launch
svr_setting = TOML.parsefile(SERVER_SETTING_FILE)
nodes = host_load_node()

planed_script = Tuple{String,Second}[]
t_long = Second(6)
t_short = Second(2)
if svr_setting["hostname"] == nodes.host.hostname
    append!(planed_script, [
        ("host/submit_to_host_buffer.jl", t_long),
        ("host/upload_to_server.jl", t_long),
        ("host/download_result.jl", t_long),
    ])
end

append!(planed_script, [
    ("server/unpack_input_file.jl", t_long),
    ("server/call_inverse.jl", t_long),
    ("server/pack_result.jl", t_long),
    ("server/update_server_status.jl", t_short)
])

if DEBUG
    @info "Planned scripts"
    @info planed_script
end

lastrun = fill(now() - Second(1) - maximum(getindex.(planed_script, 2)), length(planed_script))
fstop = joinpath(@__DIR__, hashstr(abspath(@__FILE__))*".stop")
while true
    global lastrun
    local t = now()
    if isfile(fstop)
        @info "[$t] Stopping"
        rm(fstop)
        break
    end
    for i = eachindex(planed_script)
        if t - lastrun[i] >= planed_script[i][2]
            @info "[$t] run script: $(planed_script[i][1])"
            local cmd = Cmd(`julia $(planed_script[i][1])`; dir=@__DIR__)
            # run(cmd, devnull, devnull, devnull; wait=false)
            lastrun[i] = t
        end
    end
    sleep(1)
end