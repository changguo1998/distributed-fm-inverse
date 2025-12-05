#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

nodes = host_load_node()

for s in nodes.servers
    try
        run(Cmd([
            "ssh",
            s.user*"@"*s.ip,
            """
echo 'hostname = "$(s.hostname)"' > $(replace(SERVER_SETTING_FILE, PRJ_ROOT_PATH=>s.system_root))
echo "max_event_number = $(s.max_event_number)" >> $(replace(SERVER_SETTING_FILE, PRJ_ROOT_PATH=>s.system_root))
echo "threads_per_event = $(s.threads_per_event)" >> $(replace(SERVER_SETTING_FILE, PRJ_ROOT_PATH=>s.system_root))
            """
        ]))
        @info "server $(s.hostname) updated"
    catch
        continue
    end
end