#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

if isempty(ARGS)
    println("Usage: julia update_script_file.jl <script_file>")
    exit(0)
end

nodes = host_load_node()

fpath = joinpath(PRJ_ROOT_PATH, ARGS[1])

if !isfile(fpath)
    @error "File not found: $fpath"
    exit(0)
end

for s in nodes.servers
    try
        run(Cmd(["scp", fpath, "$(s.user)@$(s.ip):$(joinpath(s.system_root, ARGS[1]))"]), devnull, devnull, devnull)
        @info "server $(s.hostname) updated"
    catch err
        @error "failed to update server $(s.hostname)"
        @error err
    end
end