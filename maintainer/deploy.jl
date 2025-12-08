#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

nodes = host_load_node()

function help()
    println("""
    Usage: julia deploy.jl hostname

    available host names:""")
    for s in nodes.servers
        println(repeat(" ", 8), s.hostname)
    end
    return nothing
end

if (length(ARGS) < 1) || ("-h" in ARGS) || ("--help" in ARGS)
    help()
    exit(0)
end
target_hostname = ARGS[1]
idxs = findall(s->s.hostname==target_hostname, nodes.servers)
if isempty(idxs)
    @info "host name not found"
    help()
    exit(0)
end
svr = nodes.servers[first(idxs)]
cmd = Cmd(["rsync", "-az", "--delete", "--delete-excluded",
    "--exclude=log/",
    "--exclude=var/",
    "--exclude=test/",
    "--exclude=julia/",
    abspath(@__DIR__, ".."),
    svr.user*"@"*svr.ip*":"*svr.system_root
    ])

run(cmd, devnull, devnull, devnull)

run(Cmd([
    "ssh",
    svr.user*"@"*svr.ip,
    """
echo 'hostname = "$(svr.hostname)"' > $(SERVER_SETTING_FILE(svr))
echo "max_event_number = $(svr.max_event_number)" >> $(SERVER_SETTING_FILE(svr))
echo "threads_per_event = $(svr.threads_per_event)" >> $(SERVER_SETTING_FILE(svr))
    """
]))
@info "Deployed to $target_hostname successfully"