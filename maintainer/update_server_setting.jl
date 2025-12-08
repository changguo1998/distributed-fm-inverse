#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

nodes = host_load_node()

for s in nodes.servers
    try
        run(Cmd([
            "ssh",
            s.user*"@"*s.ip,
            """
echo 'hostname = "$(s.hostname)"' > $(SERVER_SETTING_FILE(svr))
echo "max_event_number = $(s.max_event_number)" >> $(SERVER_SETTING_FILE(svr))
echo "threads_per_event = $(s.threads_per_event)" >> $(SERVER_SETTING_FILE(svr))
            """
        ]))
        @info "server $(s.hostname) updated"
    catch
        continue
    end
end