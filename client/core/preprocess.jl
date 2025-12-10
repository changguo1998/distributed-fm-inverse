using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using DelimitedFiles, SeisTools, TOML, Dates, JuliaSourceMechanism,
    Printf, JLD2, Statistics, LinearAlgebra

Geo = SeisTools.Geodesy

include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "multistage_lib.jl"))
include(joinpath(@__DIR__, "NLLOCshell.jl"))
include(joinpath(@__DIR__, "glibio.jl"))

function _int32_to_flag(num::Int32)
    u8num = reinterpret(reshape, UInt8, [num])
    u8str = string.(u8num, base=16, pad=2)
    return uppercase(join(u8str))
end

misfitmodules = [XCorr, Polarity]

RTSv1_setting = TOML.parsefile(joinpath(@__DIR__, "../../setting_v1.toml"))

# eventname = "202410011634"
eventname = ARGS[1]

eventpath = abspath(RTSv1_setting["focalmechanism"]["dir"], eventname)

if !isdir(eventpath)
    @error "Event $eventname not exist"
    exit(0)
end
# (_, eventname) = splitdir(abspath(eventpath))

dataroot = eventpath

if RTSv1_setting["focalmechanism"]["greenfun_type"] != "2d_compressed"
    error("Green function type not supported")
end

glibdir = abspath(RTSv1_setting["GreenFunction"]["dir"], RTSv1_setting["focalmechanism"]["greenfun_name"])
glibpath = abspath(glibdir, "glib_1.bin")
glibsetting = TOML.parsefile(joinpath(glibdir, "setting.toml"))
glib_nr = glibsetting["receiver"][1]["n"]
glib_r = range(glib_nr[1], glib_nr[3]; step=glib_nr[2])
(glib_min_dist, glib_max_dist) = extrema(abs.(glib_r .- glibsetting["srcloc"][1]))
glib_min_depth = glibsetting["receiver"][1]["d"][1]
glib_max_depth = glibsetting["receiver"][1]["d"][3]

event = TOML.parsefile(joinpath(eventpath, "event.toml"))
if event["depth"] < glib_min_depth || event["depth"] > glib_max_depth
    event["depth"] = (glib_min_depth + glib_max_depth) / 2.0
end

algorithm = deepcopy(RTSv1_setting["focalmechanism"]["algorithm"])
algorithm["searchdepth"] = event["depth"]
if algorithm["step2_min_depth"] < glib_min_depth
    algorithm["step2_min_depth"] = glib_min_depth
end
if algorithm["step2_max_depth"] > glib_max_depth
    algorithm["step2_max_depth"] = glib_max_depth
end
if algorithm["step3_min_depth"] < glib_min_depth
    algorithm["step3_min_depth"] = glib_min_depth
end
if algorithm["step3_max_depth"] > glib_max_depth
    algorithm["step3_max_depth"] = glib_max_depth
end

# unpack waveform data
sgyfile = joinpath(eventpath, "wave.sgy")
(sgyhdr, sgydata) = SeisTools.SEGY.read(sgyfile)
sacpath=joinpath(eventpath, "sac")
mkpath(sacpath)
for tr in sgydata
    trtag = _int32_to_flag(Int32(tr.hdr["sourceWaterDepth"]))
    if !(trtag in keys(event["phase"]))
        continue
    end
    sachdr = SeisTools.SAC.emptyheader()
    sachdr["delta"] = tr.hdr["dt"]*1.0e-6
    sachdr["npts"] = length(tr.data)
    sachdr["knetwk"] = String(trtag[1:2])
    sachdr["kstnm"] = String(trtag[5:8])
    sachdr["kinst"] = "UGL3C"
    if tr.hdr["groupWaterDepth"] == 0
        sachdr["kcmpnm"] = "SHN"
        sachdr["cmpaz"] = 0.0
        sachdr["cmpinc"] = 90.0
    elseif tr.hdr["groupWaterDepth"] == 1
        sachdr["kcmpnm"] = "SHE"
        sachdr["cmpaz"] = 90.0
        sachdr["cmpinc"] = 90.0
    elseif tr.hdr["groupWaterDepth"] == 2
        sachdr["kcmpnm"] = "SHZ"
        sachdr["cmpaz"] = 0.0
        sachdr["cmpinc"] = 0.0
    else
        error("unknown component")
    end
    sachdr["stla"] = event["phase"][trtag]["lat"]
    sachdr["stlo"] = event["phase"][trtag]["lon"]
    sachdr["stel"] = 0.0
    sachdr["nzyear"] = tr.hdr["yearDataRecorded"]
    sachdr["nzjday"] = tr.hdr["dayOfYear"]
    sachdr["nzhour"] = tr.hdr["hourOfDay"]
    sachdr["nzmin"] = tr.hdr["minuteOfHour"]
    sachdr["nzsec"] = tr.hdr["secondOfMinute"]
    sachdr["nzmsec"] = 0
    sachdr["b"] = 0.0
    sacfile = joinpath(sacpath, SeisTools.SAC.standardname(sachdr))
    if sachdr["kcmpnm"] == "SHZ"
        SeisTools.SAC.write(sacfile, sachdr, -tr.data)
    else
        SeisTools.SAC.write(sacfile, sachdr, tr.data)
    end
end

# stations = buildstationconfiguration(dataroot, event)

# filter!(stations) do s
#     return (s["base_distance"] > glib_min_dist) && (s["base_distance"] < glib_max_dist)
# end

tstations = buildstationconfiguration(dataroot, event)
filter!(tstations) do s
    return (s["base_distance"] > glib_min_dist) && (s["base_distance"] < glib_max_dist)
end
stationtag(s) = s["network"]*"."*s["station"]
stations = let
    station_list = unique(sort(map(stationtag, tstations)))
    w = randn(length(station_list))
    r = sortperm(w)
    if length(station_list) > 6
        station_selected = station_list[r[1:6]]
    else
        station_selected = station_list
    end
    filter(s->stationtag(s) in station_selected, tstations)
end

TIME_DEFAULT = DateTime(2000)

# = station infomation

for s in stations
    # println(s["station"])
    sacf = SeisTools.SAC.read(normpath(dataroot, "sac", s["meta_file"]))
    # - trim window
    reft = SeisTools.SAC.DateTime(sacf.hdr)
    bt = event["origintime"] - Minute(1)
    et = bt + Minute(5)
    # et = reft + msecond(sacf.hdr["delta"] * sacf.hdr["npts"])
    s["base_trim"] = [bt, et]
    stag = join([s["network"], "00", s["station"]])
    if stag in keys(event["phase"])
        if "P" in keys(event["phase"][stag])
            push!(s["phases"], Dict("type"=>"P", "at"=>event["phase"][stag]["P"]))
        end

        if "S" in keys(event["phase"][stag])
            push!(s["phases"], Dict("type"=>"S", "at"=>event["phase"][stag]["S"]))
        end
    else
        push!(s["phases"], Dict(
            "type" => "P",
            "at" => TIME_DEFAULT
        ))
    end

    if isnan(sacf.hdr["scale"])
        s["meta_scale"] = 1.0
    else
        s["meta_scale"] = sacf.hdr["scale"]
    end
    # - Green function setting
    # general options
    s["green_modeltype"] = "2D_COMPRESSED"
    s["green_model"] = RTSv1_setting["focalmechanism"]["greenfun_name"]
    s["green_tsource"] = glibsetting["risetime"]
    s["green_dt"] = glibsetting["dt"]
    s["green_modelpath"] = glibpath
    (_, _, tp, ts, _) = JuliaSourceMechanism.Green.cglib_readlocation(glibpath, -s["base_distance"], 0.0, event["depth"])
    # * phase infomation
    for p in s["phases"]
        # p["tt"] = round(p["at"] - event["origintime"], Millisecond).value * 1e-3
        if p["type"] == "P"
            defaultband = [0.05, 0.2]
            p["tt"] = tp
        else
            defaultband = [0.05, 0.1]
            p["tt"] = ts
        end
        p["xcorr_order"] = 4
        p["xcorr_band"] = defaultband
        p["xcorr_maxlag"] = 1.0 / p["xcorr_band"][2]
        if p["type"] == "P"
            p["xcorr_trim"] = [-2.0, 3.0] ./ p["xcorr_band"][2]
            p["polarity_obs"] = 0.0
            p["polarity_trim"] = [0.0, s["green_tsource"]]
        else
            p["xcorr_trim"] = [-4.0, 6.0] ./ p["xcorr_band"][2]
            p["polarity_obs"] = NaN
            p["polarity_trim"] = [NaN, NaN]
        end
        local dt = s["meta_dt"]
        tl = p["xcorr_trim"][2] - p["xcorr_trim"][1]
        while (dt + 1e-3) * 200 < tl
            dt += 1e-3
        end
        p["xcorr_dt"] = dt
    end
end

# = construct env data

env = Dict("algorithm" => algorithm,
        "event" => event,
        "stations" => stations,
        "dataroot" => dataroot)

JuliaSourceMechanism.loaddata!(env)
for s in env["stations"]
    ot = Millisecond(env["event"]["origintime"] - s["base_begintime"]).value * 1e-3
    (meta, _) = JuliaSourceMechanism.Green.scangreenfile(normpath(env["dataroot"],
                                                                "greenfun",
                                                                @sprintf("%s-%.4f", s["green_model"],
                                                                        env["algorithm"]["searchdepth"]),
                                                                s["network"] * "." * s["station"] * "." *
                                                                s["component"] * ".gf"))
    for p in s["phases"]
        if p["type"] == "P"
            p["tt"] = meta["tp"] + ot
        elseif p["type"] == "S"
            p["tt"] = meta["ts"] + ot
        end
    end
end
status = Dict{String,Any}()
status["saveplotdata"] = true
status["saveplotdatato"] = abspath(".tmpplot.mat")

# = copy data and delete extra station

tenv = Dict{String,Any}()
tenv["algorithm"] = deepcopy(env["algorithm"])
tenv["event"] = deepcopy(env["event"])
tenv["dataroot"] = deepcopy(env["dataroot"])
scpflag = falses(length(env["stations"]))
cpcheck = falses(length(env["stations"]))
while !all(scpflag)
    i = findfirst(!, scpflag)
    js = findall(s->(s["network"] == env["stations"][i]["network"]) &&
        (s["station"] == env["stations"][i]["station"]), env["stations"])
    w = 0.0
    for j in js
        for p in env["stations"][j]["phases"]
            w += ("xcorr_weight" in keys(p)) ? p["xcorr_weight"] : 1.0
        end
    end
    if !iszero(w)
        for j in js
            cpcheck[j] = true
        end
    end
    for j in js
        scpflag[j] = true
    end
end
tenv["stations"] = Vector{Dict{String,Any}}(undef, count(cpcheck))
cpid = 1
for i = eachindex(env["stations"])
    global cpid
    if cpcheck[i]
        tenv["stations"][cpid] = deepcopy(env["stations"][i])
        cpid += 1
    end
end

jldsave(joinpath(dataroot, "auto.jld2"); env=tenv, status=status)
