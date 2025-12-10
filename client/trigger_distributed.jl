using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using TOML, Dates, Printf

setting = TOML.parsefile(joinpath(@__DIR__, "../setting_v1.toml"))

stop_sign_file = joinpath(@__DIR__, "stop_trigger.flag")
inversion_failed_sign_filename = "inversion_failed.flag"
plot_end_flag = "plot_end.flag"
plot_failed_flag = "plot_failed.flag"
const FLAG_HOST_PREPROCESS_BEGIN = "dfmi_host_preprocess_begin.flag"
const FLAG_HOST_PREPROCESS_END = "dfmi_host_preprocess_end.flag"

while true
    stop_sign = isfile(stop_sign_file)
    if stop_sign
        @info "Stop trigger detected. Exiting script."
        rm(stop_sign_file)
        break
    end

    events = readdir(setting["focalmechanism"]["dir"])

    events_not_inversed = filter(events) do e
        epath = joinpath(setting["focalmechanism"]["dir"], e)
        erpath = joinpath(epath, "result")
        failed_flag = isfile(joinpath(epath, inversion_failed_sign_filename))
        return !isdir(erpath) &&
            !failed_flag &&
            !isfile(epath, FLAG_HOST_PREPROCESS_BEGIN)
    end
    if !isempty(events_not_inversed)
        event_scheduled = events_not_inversed[end]

        @info "Inverse event $event_scheduled"
        try
            run(Cmd(
                Cmd(["julia", "inverse_caller_distributed.jl", event_scheduled])
                ; dir=@__DIR__)
            )
        catch err
            touch(joinpath(setting["focalmechanism"]["dir"],
                event_scheduled,
                inversion_failed_sign_filename))
        end
        @info "Done $event_scheduled"
    end

    events_not_plot = filter(events) do e
        epath = joinpath(setting["focalmechanism"]["dir"], e)
        erpath = joinpath(epath, "result")
        if !isdir(erpath)
            return false
        end
        if isfile(joinpath(epath, inversion_failed_sign_filename))
            return false
        end
        if isfile(joinpath(epath, plot_failed_flag))
            return false
        end
        if isfile(joinpath(epath, plot_end_flag))
            return false
        end
        n_jld2 = length(filter(endswith(".jld2"), readdir(erpath)))
        n_png = length(filter(endswith(".png"), readdir(erpath)))
        return (n_jld2 == 3) && iszero(n_png)
    end
    if !isempty(events_not_plot)
        event_sc = events_not_plot[end]
        @info "Plot event: $event_sc"
        try
            run(Cmd(
                Cmd(["julia", "plot_caller.jl", event_sc])
                ; dir=@__DIR__)
            )
            touch(joinpath(setting["focalmechanism"]["dir"],
                event_sc,
                plot_end_flag))
        catch
            touch(joinpath(setting["focalmechanism"]["dir"],
                event_sc,
                plot_failed_flag))
        end
    end
    # run(Cmd(Cmd(["julia", "collect_result.jl"]); dir=@__DIR__))
    sleep(setting["focalmechanism"]["trigger_time_interval"])
end
