#!/usr/bin/env julia

include(joinpath(@__DIR__, "../lib.jl"))

svr_setting = TOML.parsefile(SERVER_SETTING_FILE)
nodes = host_load_node()

if svr_setting["hostname"] == nodes.host.hostname
end