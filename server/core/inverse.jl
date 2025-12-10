using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using DelimitedFiles, SeisTools, TOML, Dates, JuliaSourceMechanism,
    Printf, JLD2, Statistics, LinearAlgebra

Geo = SeisTools.Geodesy

include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "multistage_lib.jl"))
include(joinpath(@__DIR__, "NLLOCshell.jl"))
include(joinpath(@__DIR__, "glibio.jl"))

# usage:
#     julia inverse.jl eventname

eventpath = abspath(ARGS[1])

(env, status) = let
    datafilepath = joinpath(eventpath, "auto.jld2")
    if !isfile(datafilepath)
        @error "File $datafilepath not exist"
        exit(0)
    end
    # load data
    t = load(datafilepath)
    (t["env"], t["status"])
end

# @info "Copy data"
rawenv = Setting()
rawenv["dataroot"] = eventpath
dataroot=rawenv["dataroot"]
rawenv["algorithm"] = env["algorithm"]
rawenv["event"] = env["event"]
rawenv["stations"] = Setting[]
for s in env["stations"]
    t = deepcopy(s)
    push!(rawenv["stations"], t)
end

misfits = Module[]
for m in rawenv["algorithm"]["misfit"], f in [XCorr, Polarity, PSR, DTW, AbsShift, RelShift]
    if m in f.tags
        push!(misfits, f)
    end
end

if length(rawenv["stations"]) < rawenv["algorithm"]["minimum_stations"]
    error("No enough stations")
end

# @info "Run"
JuliaSourceMechanism.calcgreen!(rawenv)

(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(rawenv, misfits)

nenv = deepcopy(rawenv)
hlist = [nenv["algorithm"]["searchdepth"]]
misfitlist = [minval]
@write_result "result_stage1"

(mech, minval, minval_xcorr, _minval_pol) = let
    local _tfms = zeros(3, 0)
    local _mech = zeros(3)
    local _tmech = zeros(3)
    _mech .= mech
    local _minval = 0.0
    local _minval_xcorr = 0.0
    local _minval_pol = 0.0
    local Nsample = 3
    local Ncount = 0
    while true
        _rmech = randn(3)
        _rmech[1] = min(3.0, max(-3.0, _rmech[1]))
        _rmech[2] = min(3.0, max(-3.0, _rmech[2]))
        _rmech[3] = min(3.0, max(-3.0, _rmech[3]))
        _tmech[1] = mod(_mech[1]+_rmech[1], 360.0)
        _tmech[2] = max(0.0, min(90.0, _mech[2]+_rmech[2]))
        _tmech[3] = min(90.0, max(-90.0, _mech[3] + _rmech[3]))
        update_filter_band!(nenv, dc2ts(_tmech), 0.2:0.01:0.25, "P", -2.0, 3.0)
        update_filter_band!(nenv, dc2ts(_tmech), 0.1:0.01:0.15, "S", -4.0, 6.0)
        (_mech, _minval, _minval_xcorr, _minval_pol) = inverse_focalmech!(nenv, misfits)
        _tfms = hcat(_tfms, _mech)
        if size(_tfms, 2) < Nsample
            # println(mech)
            continue
        end
        cstrike = cosd.(_tfms[1,end-Nsample+1:end])
        sstrike = sind.(_tfms[1,end-Nsample+1:end])
        std1 = sqrt(std(cstrike)^2+std(sstrike)^2)
        std2 = std(_tfms[2,end-Nsample+1:end])
        std3 = std(_tfms[3,end-Nsample+1:end])
        # println([mech; std1; std2; std3])
        if std1 < sind(5) && max(std2, std3) < 2.5
            break
        end
        Ncount += 1
        if Ncount > nenv["algorithm"]["frequency_test_maximum_iteration"]
            break
        end
    end
    (_mech, _minval, _minval_xcorr, _minval_pol)
end

hmin = nenv["algorithm"]["step2_min_depth"]
dh = nenv["algorithm"]["step2_d_depth"]
hmax = nenv["algorithm"]["step2_max_depth"]
r = nenv["algorithm"]["step2_depth_radius"]

h = round((nenv["algorithm"]["searchdepth"] - hmin)/dh) * dh + hmin
hlist = Float64[]
misfitlist = Float64[]
for _ = 1:10
    global h
    local tl = filter(_h->!(_h in hlist), max(hmin, h-r):dh:min(hmax, h+r))
    local val = inverse_depth(tl, nenv, misfits)
    append!(hlist, tl)
    append!(misfitlist, val)
    local h2 = hlist[argmin(misfitlist)]
    if abs(h2 - h) < nenv["algorithm"]["step2_stop_iterate_when_depth_change_less_than"]
        h = h2
        break
    end
    h = h2
end

nenv["algorithm"]["searchdepth"] = h
(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

@write_result "result_stage2"

thres_list = 0.0:0.1:1.0
select_result = map(thres_list) do thres
    local test_env = reselect_channel(nenv, dc2ts(mech), thres)
    return length(test_env["stations"])
end
thresIdx = findlast(select_result) do n
    if length(nenv["stations"]) < 3
        return true
    elseif length(nenv["stations"]) < 7
        return n > 3
    else
        return n >= 6
    end
end
xcorr_threshold = thres_list[thresIdx]
nenv = reselect_channel(nenv, dc2ts(mech), xcorr_threshold)
nenv["algorithm"]["xcorr_threshold"] = xcorr_threshold
# (mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

hmin = nenv["algorithm"]["step3_min_depth"]
dh = nenv["algorithm"]["step3_d_depth"]
hmax = nenv["algorithm"]["step3_max_depth"]
r = nenv["algorithm"]["step3_depth_radius"]

h = round((nenv["algorithm"]["searchdepth"] - hmin)/dh) * dh + hmin
hlist = Float64[]
misfitlist = Float64[]
for _ = 1:10
    global h
    local tl = filter(_h->!(_h in hlist), max(hmin, h-r):dh:min(hmax, h+r))
    local val = inverse_depth(tl, nenv, misfits)
    append!(hlist, tl)
    append!(misfitlist, val)
    local h2 = hlist[argmin(misfitlist)]
    if abs(h2 - h) < nenv["algorithm"]["step3_stop_iterate_when_depth_change_less_than"]
        h = h2
        break
    else
        h = h2
    end
end

nenv["algorithm"]["searchdepth"] = h
(mech, minval, minval_xcorr, minval_pol) = inverse_focalmech!(nenv, misfits)

@write_result "result_stage3"
