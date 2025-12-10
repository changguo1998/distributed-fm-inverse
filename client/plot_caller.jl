using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using TOML

setting = TOML.parsefile(joinpath(@__DIR__, "../setting_v1.toml"))

eventname = ARGS[1]
eventpath = joinpath(setting["focalmechanism"]["dir"], eventname)

if !isdir(eventpath)
    @error "Event $eventname not exist"
    exit(0)
end
# exit(0)
run(Cmd(Cmd([
    "julia",
    "inv_script/plot.jl",
    eventname
]); dir=@__DIR__))
