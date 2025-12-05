#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "../lib.jl"))

nodes = host_load_node()

function help()
    println("""
    Usage: julia deploy.jl hostname

    available host names:
""")
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
    help()
    exit(0)
end
svr = hosts.servers[first(idxs)]\
cmd = Cmd(`rsync -avz --delete $(abspath(@__DIR__, "..")) $(svr.user)@$(svr.ip):$(svr.system_root)`)

run(cmd, devnull, devnull, devnull)
@info "Succeed to deploy to $target_hostname"