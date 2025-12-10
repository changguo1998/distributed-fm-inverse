#!/usr/bin/env julia
include(joinpath(@__DIR__, "lib.jl"))

const FAST_FLUSH = 1
const SLOW_FLUSH = 2

function test_loading_time()
    @info "Test loading time"
    t1 = now()
    t2 = t1
    try
        run(Cmd(`bash dfmi.sh loading_time_test.jl`; dir=@__DIR__), devnull, devnull, devnull)
    finally
        t2 = now()
    end
    return max(ceil(t2 - t1, Millisecond), Millisecond(1000))
end

# reset
cmd = Cmd(`bash dfmi.sh reset.jl`; dir=@__DIR__)
if DEBUG
    @info cmd
end
run(cmd, devnull, devnull, devnull)
run(`bash dfmi.sh server/update_server_status.jl`, devnull, devnull, devnull)

# launch
svr_setting = TOML.parsefile(SERVER_SETTING_FILE())
nodes = host_load_node()
loading_time = let
    n = 5
    Millisecond(round(Int,
        sum(map(i->test_loading_time(), 1:n)).value/n
    ))
end

planed_script = Tuple{String,Second,Int}[]
interval = max(ceil(loading_time, Second), Second(10))
@info "run time interval: $interval"
if svr_setting["hostname"] == nodes.host.hostname
    append!(planed_script, [
        ("host/submit_to_host_buffer.jl", interval, SLOW_FLUSH),
        ("host/upload_to_server.jl", interval, SLOW_FLUSH),
        ("host/download_result.jl", interval, SLOW_FLUSH),
    ])
    idxs = findfirst(s->s.hostname == nodes.host.hostname, nodes.servers)
    if !isnothing(idxs)
        if nodes.servers[idxs].priority >= 0
            append!(planed_script, [
                ("server/unpack_input_file.jl", interval, SLOW_FLUSH),
                ("server/call_inverse.jl", interval, SLOW_FLUSH),
                ("server/pack_result.jl", interval, SLOW_FLUSH),
                ("server/update_server_status.jl", interval, FAST_FLUSH),
                ("server/cleanup_downloaded_result.jl", interval, FAST_FLUSH)
            ])
        end
    end
else
    append!(planed_script, [
        ("server/unpack_input_file.jl", interval, SLOW_FLUSH),
        ("server/call_inverse.jl", interval, SLOW_FLUSH),
        ("server/pack_result.jl", interval, SLOW_FLUSH),
        ("server/update_server_status.jl", interval, FAST_FLUSH),
        ("server/cleanup_downloaded_result.jl", interval, FAST_FLUSH)
    ])
end

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
    slow_flush_script_is_called = false
    for i = eachindex(planed_script)
        t = now()
        script = planed_script[i]
        if ((t - lastrun[i] >= script[2]) && script[3]==SLOW_FLUSH) || (script[3]==FAST_FLUSH)
            if DEBUG
                @info "[$t] run script: $(script[1])"
            end
            local cmd = Cmd(`bash dfmi.sh $(script[1])`; dir=@__DIR__)
            run(cmd; wait=false)
            lastrun[i] = t
            if script[3] == SLOW_FLUSH
                slow_flush_script_is_called = true
            end
        end
    end
    if !slow_flush_script_is_called
        sleep(loading_time)
    end
end