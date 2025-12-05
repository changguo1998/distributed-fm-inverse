#!/usr/bin/env julia

include(joinpath(@__DIR__, "lib.jl"))

mkpath(BUFFER_HOST_UPLOAD)
mkpath(BUFFER_HOST_RESULT)
mkpath(BUFFER_SERVER_INPUT)
mkpath(BUFFER_SERVER_RUN)
mkpath(BUFFER_SERVER_RESULT)
mkpath(joinpath(PRJ_ROOT_PATH, "log"))