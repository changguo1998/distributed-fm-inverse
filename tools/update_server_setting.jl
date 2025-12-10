#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

nodes = host_load_node()

for s in nodes.servers
    try
        run(Cmd(["scp", NODE_LIST_FILE(), "$(s.user)@$(s.ip):$(NODE_LIST_FILE(s))"]), devnull, devnull, devnull)
        run(Cmd([
            "ssh",
            s.user*"@"*s.ip,
            """
echo 'hostname = "$(s.hostname)"' > $(SERVER_SETTING_FILE(s))
echo "max_event_number = $(s.max_event_number)" >> $(SERVER_SETTING_FILE(s))
echo "threads_per_event = $(s.threads_per_event)" >> $(SERVER_SETTING_FILE(s))
            """
        ]))
        @info "server $(s.hostname) updated"
    catch err
        @error "failed to update server $(s.hostname)"
        @error err
        continue
    end
end