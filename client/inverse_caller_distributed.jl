using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
using TOML
const FLAG_HOST_PREPROCESS_BEGIN = "dfmi_host_preprocess_begin.flag"
const FLAG_HOST_PREPROCESS_END = "dfmi_host_preprocess_end.flag"
setting = TOML.parsefile(joinpath(@__DIR__, "../setting_v1.toml"))

eventname = ARGS[1]
eventpath = joinpath(setting["focalmechanism"]["dir"], eventname)

if !isdir(eventpath)
    @error "Event $eventname not exist"
    exit(0)
end
touch(joinpath(eventpath, FLAG_HOST_PREPROCESS_BEGIN))
run(Cmd(Cmd(["julia", "inv_script/preprocess_distributed.jl", eventname]); dir=@__DIR__))
run(Cmd(Cmd([
    "julia",
    "inv_script/inverse_distributed.jl",
    eventname
    ]); dir=@__DIR__))
touch(joinpath(eventpath, FLAG_HOST_PREPROCESS_END))
