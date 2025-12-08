#!/usr/bin/env julia
include(joinpath(@__DIR__, "lib.jl"))

fstop = joinpath(@__DIR__, hashstr(abspath(@__DIR__, "launch.jl"))*".stop")

touch(fstop)

while isfile(fstop)
    sleep(0.5)
end

@info "monitor stopped"

@info "backup started"
bkfile = "bak_"*Dates.format(now(), "yyyymmddHHMMSS")*".tar.gz"
cmd = Cmd(`tar czf $bkfile log var`; dir=@__DIR__)
if DEBUG
    @info cmd
end
run(cmd, devnull, devnull, devnull)
@info "backup data to $bkfile"
