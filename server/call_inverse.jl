#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))
const LOG_SETTING = (log=LOG_SERVER_INVERSE(), lock=LOCK_SERVER_INVERSE_LOG)

function main()
    buffer_run_list = readdir(BUFFER_SERVER_RUN())
    if isempty(buffer_run_list)
        return nothing
    end
    wait_for_submit = filter(buffer_run_list) do d
        edir = joinpath(BUFFER_SERVER_RUN(), d)
        return isfile(joinpath(edir, FLAG_SERVER_UNPACKED)) && (!isfile(joinpath(edir, FLAG_SERVER_INVERSE_BEGIN)))
    end
    if isempty(wait_for_submit)
        return nothing
    end
    evt = first(wait_for_submit)
    log_info("call inverse for event $evt")
    touch(joinpath(BUFFER_SERVER_RUN(), evt, FLAG_SERVER_INVERSE_BEGIN))
    svrsetting = TOML.parsefile(SERVER_SETTING_FILE())
    callfile = DRY_RUN ? "inverse_dry.jl" : "inverse.jl"
    tenv = deepcopy(ENV)
    tenv["JULIA_DEPOT_PATH"] = abspath(PRJ_ROOT_PATH, "julia")
    cmd = Cmd(Cmd([joinpath(PRJ_ROOT_PATH, "julia", "bin", "julia"), "-t", string(svrsetting["threads_per_event"]),
        joinpath(@__DIR__, "core", callfile), joinpath(BUFFER_SERVER_RUN(), evt)]); env = tenv)

    try
        run(cmd)
        touch(joinpath(BUFFER_SERVER_RUN(), evt, FLAG_SERVER_INVERSE_END))
        log_info("inversion for event $evt succeeded")
    catch err
        touch(joinpath(BUFFER_SERVER_RUN(), evt, FLAG_SERVER_INVERSE_FAILED))
        log_err("inversion for event $evt failed")
        error(err)
    end
end

try
    main()
catch err
    log_err(string(err))
end