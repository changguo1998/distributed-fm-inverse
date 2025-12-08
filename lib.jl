# ENV["JULIA_DEPOT_PATH"] = abspath(@__DIR__, "julia/")
using Pkg
Pkg.activate(@__DIR__)

using TOML, Dates, SHA

const PRJ_ROOT_PATH = abspath(@__DIR__)
const DRY_RUN = true
const DEBUG = true

# buffer dir
const BUFFER_HOST_UPLOAD = abspath(@__DIR__, "var/host_upload")
const BUFFER_HOST_RESULT = abspath(@__DIR__, "var/host_result")
const BUFFER_SERVER_INPUT = abspath(@__DIR__, "var/svr_input")
const BUFFER_SERVER_RUN = abspath(@__DIR__, "var/svr_run")
const BUFFER_SERVER_RESULT = abspath(@__DIR__, "var/svr_result")

# flag file
# client event dir
const FLAG_HOST_PREPROCESS_BEGIN = "dfmi_host_preprocess_begin.flag"
const FLAG_HOST_PREPROCESS_END = "dfmi_host_preprocess_end.flag"
const FLAG_HOST_QUEUE_BEGIN = "dfmi_host_queue_begin.flag"
const FLAG_HOST_QUEUE_END = "dfmi_host_queue_end.flag"
const FLAG_HOST_INVERSION_FINISHED = "dfmi_host_inversion_finished.flag"

# server event dir
const FLAG_SERVER_UNPACKED = "dfmi_svr_unpacked.flag"
const FLAG_SERVER_INVERSE_BEGIN = "dfmi_svr_inverse_begin.flag"
const FLAG_SERVER_INVERSE_FAILED = "dfmi_svr_inversion_failed.flag"
const FLAG_SERVER_INVERSE_END = "dfmi_svr_inverse_end.flag"
const FLAG_SERVER_PACK_RESULT_BEGIN = "dfmi_svr_pack_result_begin.flag"
const FLAG_SERVER_PACK_RESULT_END = "dfmi_svr_pack_result_end.flag"

# server var
const FLAG_SERVER_UPLOADED = abspath(PRJ_ROOT_PATH, "var/dfmi_svr_uploaded.flag")

# log file
const LOG_HOST_QUEUE = abspath(PRJ_ROOT_PATH, "log/host_submit.log")
const LOG_HOST_UPLOAD = abspath(PRJ_ROOT_PATH, "log/host_upload.log")
const LOG_HOST_DOWNLOAD = abspath(PRJ_ROOT_PATH, "log/host_download.log")
const LOG_SERVER_UNPACK = abspath(PRJ_ROOT_PATH, "log/server_unpack.log")
const LOG_SERVER_INVERSE = abspath(PRJ_ROOT_PATH, "log/server_inverse.log")
const LOG_SERVER_RESULT = abspath(PRJ_ROOT_PATH, "log/server_result.log")

# status file
const STATUS_SERVER = abspath(PRJ_ROOT_PATH, "var/svr_status.toml")
const STATUS_HOST = abspath(PRJ_ROOT_PATH, "var/host_status.toml")
const STATUS_QUEUE = abspath(PRJ_ROOT_PATH, "var/queue_status.toml")

# lock file
const LOCK_HOST_QUEUE_STATUS_FILE = abspath(PRJ_ROOT_PATH, "var/queue_status.lock")

const LOCK_HOST_QUEUE_LOG = abspath(PRJ_ROOT_PATH, "var/host_queue_log.lock")
const LOCK_HOST_UPLOAD_LOG = abspath(PRJ_ROOT_PATH, "var/host_upload_log.lock")
const LOCK_HOST_DOWNLOAD_LOG = abspath(PRJ_ROOT_PATH, "var/host_download_log.lock")
const LOCK_SERVER_UNPACK_LOG = abspath(PRJ_ROOT_PATH, "var/server_unpack_log.lock")
const LOCK_SERVER_INVERSE_LOG = abspath(PRJ_ROOT_PATH, "var/server_inverse_log")
const LOCK_SERVER_RESULT_LOG = abspath(PRJ_ROOT_PATH, "var/server_result_log.lock")

const NODE_LIST_FILE = abspath(PRJ_ROOT_PATH, "config/node-list-test-host-server.toml")
const SERVER_SETTING_FILE = abspath(PRJ_ROOT_PATH, "config/svr.toml")

function hashstr(s::String)
    h = sha256(s)
    return uppercase(join(bytes2hex.(h[1:8])))
end

function get_lock(f::String)
    n = 1
    while isfile(f)
        sleep(0.1 * n)
        n = rand(1:5)
    end
    return nothing
end

function release_lock(f::String)
    if isfile(f)
        rm(f)
    end
    return nothing
end

function get_single_process_lock(f::String)
    lockfile = joinpath(PRJ_ROOT_PATH, "var", hashstr(f)*".lock")
    if isfile(lockfile)
        exit(0)
    end
    touch(lockfile)
    return nothing
end

function release_single_process_lock(f::String)
    lockfile = joinpath(PRJ_ROOT_PATH, "var", hashstr(f)*".lock")
    if isfile(lockfile)
        rm(lockfile; force=true)
    end
end

function log_msg(prefix::String, msg...)
    global LOG_SETTING

    b = join(msg)
    if isempty(b)
        return nothing
    end
    lines = split(b, '\n')

    get_lock(LOG_SETTING.lock)
    t = now()
    tstr = Dates.format(t, "yyyy-mm-dd HH:MM:SS.sss")
    headline = "[$tstr|$prefix] "
    L = length(headline)
    open(LOG_SETTING.log, "a") do io
        println(io, headline, lines[1])
        for line in lines[2:end]
            println(io, " "^L, line)
        end
    end
    release_lock(LOG_SETTING.lock)
    return nothing
end

log_info(msg...) = log_msg(" INFO", msg...)
log_warn(msg...) = log_msg(" WARN", msg...)
log_err(msg...) = log_msg("ERROR", msg...)


struct InvServer <: Any
    hostname::String
    ip::String
    user::String
    system_root::String
    max_event_number::Int
    threads_per_event::Int
    priority::Int
end

function InvServer(
    hostname::AbstractString,
    ip::AbstractString,
    user::AbstractString,
    system_root::AbstractString,
    max_event_number::Integer,
    threads_per_event::Integer,
    priority::Integer)

    return InvServer(
        String(hostname),
        String(ip),
        String(user),
        String(system_root),
        Int(max_event_number),
        Int(threads_per_event),
        Int(priority)
    )
end

struct InvHost <: Any
    hostname::String
    ip::String
    user::String
    system_root::String
    monitor_directory::Vector{String}
end

function InvHost(
    hostname::AbstractString,
    ip::AbstractString,
    user::AbstractString,
    system_root::AbstractString,
    monitor_directory::AbstractVector{<:AbstractString})

    return InvHost(
        String(hostname),
        String(ip),
        String(user),
        String(system_root),
        map(String, monitor_directory)
    )
end

function host_load_node()
    t = TOML.parsefile(NODE_LIST_FILE)
    svrs = map(t["server"]) do s
        return InvServer(
            s["hostname"],
            s["ip"],
            s["user"],
            s["system_root"],
            s["max_event_number"],
            s["threads_per_event"],
            s["priority"]
        )
    end
    h = InvHost(
        t["host"]["hostname"],
        t["host"]["ip"],
        t["host"]["user"],
        t["host"]["system_root"],
        t["host"]["monitor_directory"]
    )
    return (host=h, servers=svrs)
end

function get_server_loading(svr::InvServer, host::InvHost)
    # flush
    cmd_flush = Cmd([
        "ssh",
        svr.user*"@"*svr.ip,
        "cd $(svr.system_root); bash dfmi.sh server/update_server_status.jl"
    ])
    try
        @info cmd_flush
        run(cmd_flush)
    catch
        return nothing
    end

    # download command
    if svr.hostname == host.hostname
        cmd = Cmd(["cat", STATUS_SERVER])
    else
        server_status_file = replace(STATUS_SERVER, PRJ_ROOT_PATH=>svr.system_root)
        cmd = Cmd(`ssh $(svr.user*"@"*svr.ip) "cat $(server_status_file)"`)
    end
    try
        s = read(cmd, String)
        return TOML.parse(s)
    catch
        return nothing
    end
    return nothing
end
