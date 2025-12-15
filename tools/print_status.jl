#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

using ArgumentProcessor

function main()

    ArgumentProcessor.THROW_ERROR_AFTER_HELP = true

    addflag!("all";         help="print all information")
    addflag!("queue";       abbr="q",  help="events submitted to queue")
    addflag!("config";      abbr="c",  help="server configuration")
    addflag!("loading";     abbr="l",  help="server loading status")
    addflag!("log_queue";   abbr="Lq", outername="log-queue",    help="log file of submit_to_host_buffer.jl")
    addflag!("log_up";      abbr="Lu", outername="log-upload",   help="log file of upload_to_server.jl")
    addflag!("log_down";    abbr="Ld", outername="log-download", help="log file of download_result.jl")
    addflag!("log_unpack";  abbr="Lp", outername="log-unpack",   help="log file of unpack_input_file.jl")
    addflag!("log_inverse"; abbr="Li", outername="log-inverse",  help="log file of call_inverse.jl")
    addflag!("log_result";  abbr="Lr", outername="log-result",   help="log file of pack_result.jl")

    addopt!("server";    abbr="s", fmt=" %s", default=" ", help="specify server host name while printing remote log files")
    addopt!("log_lines"; abbr="n", outername="log-lines", fmt=" %d", default=" 5", help="print last number of lines while printing log file")

    input = ArgumentProcessor.parse(ARGS)
    nodes = host_load_node()
    server_list = map(s->s.hostname, nodes.servers)
    line_sep = 0
    skip_n_line = 1
    if input.config || input.all
        if isempty(nodes.servers)
            print(repeat("\n", line_sep))
            @warn "Server list is empty"
        else
            info = fill("", length(server_list)+1, 5)
            info[1, :] = ["Server", "Priority", "Max Loading", "Address", "Root path"]
            for i in eachindex(nodes.servers)
                s = nodes.servers[i]
                info[i+1, 1] = s.hostname
                info[i+1, 2] = string(s.priority)
                info[i+1, 3] = string(s.max_event_number) * " x " * string(s.threads_per_event)
                info[i+1, 4] = s.user * "@" * s.ip
                info[i+1, 5] = s.system_root
            end
            print(repeat("\n", line_sep))
            println("Setting:")
            print_table(info; aligncenter=[2, 3])
        end
        line_sep += skip_n_line - line_sep
    end

    if input.queue || input.all
        if isfile(STATUS_QUEUE)
            qinfo = TOML.parsefile(STATUS_QUEUE)
            id_list = setdiff(collect(keys(qinfo)), ["update_time"])
            filter!(t->isfile(joinpath(BUFFER_HOST_UPLOAD, t*"_input.tar.gz")), id_list)
            if isempty(id_list)
                print(repeat("\n", line_sep))
                printstyled("Event queue is empty\n"; color=:light_green)
            else
                info = vcat(["Tag" "Path"], [isone(j) ? id_list[i] : qinfo[id_list[i]] for i = eachindex(id_list), j = 1:2])
                print(repeat("\n", line_sep))
                println("Queue:")
                print_table(string.(info))
            end
        else
            print(repeat("\n", line_sep))
            @warn "Queue status file not found."
        end
        line_sep += skip_n_line - line_sep
    end

    if input.loading || input.all
        if isempty(server_list)
            print(repeat("\n", line_sep))
            @info "No available server."
        else
            info = fill("", length(server_list)+1, 4)
            info[1, :] = ["Server", "Priority", "Loading", "Update Time"]
            for i in eachindex(nodes.servers)
                s = nodes.servers[i]
                loading_status = get_server_loading(s, nodes.host)
                info[i+1, 1] = s.hostname
                info[i+1, 2] = string(s.priority)
                if isnothing(loading_status)
                    info[i+1, 3] = "ERROR"
                    info[i+1, 4] = ""
                else
                    info[i+1, 3] = string(length(setdiff(loading_status["inverse"], loading_status["result"])))*"/"*string(s.max_event_number)
                    info[i+1, 4] = string(loading_status["update_time"])
                end
            end
            print(repeat("\n", line_sep))
            println("Loading Status:")
            print_table(info; aligncenter=[2, 3])
        end
        line_sep += skip_n_line - line_sep
    end

    if input.log_queue || input.all
        print(repeat("\n", line_sep))
        print_log(LOG_HOST_QUEUE, input.log_lines)
        line_sep += skip_n_line - line_sep
    end
    if input.log_up || input.all
        print(repeat("\n", line_sep))
        print_log(LOG_HOST_UPLOAD, input.log_lines)
        line_sep += skip_n_line - line_sep
    end
    if input.log_down || input.all
        print(repeat("\n", line_sep))
        print_log(LOG_HOST_DOWNLOAD, input.log_lines)
        line_sep += skip_n_line - line_sep
    end

    remote_info_list = String[]
    @debug "remote_info_list: $(remote_info_list)"
    if input.all
        append!(remote_info_list, server_list)
    end
    @debug "remote_info_list: $(remote_info_list)"
    if !isnothing(input.server) && !isempty(input.server)
        push!(remote_info_list, input.server)
    end
    unique!(remote_info_list)

    if isempty(remote_info_list)
        return nothing
    end

    @debug "remote_info_list: $(remote_info_list)"

    for svrhostname = remote_info_list
        i = findfirst(==(svrhostname), server_list)
        if isnothing(i)
            @error "Hostname $(svrhostname) not found in server list"
            continue
        end
        s = nodes.servers[i]
        if input.log_unpack || input.all
            print(repeat("\n", line_sep))
            print_remote_log(s, LOG_SERVER_UNPACK, input.log_lines)
            line_sep += skip_n_line - line_sep
        end
        if input.log_inverse || input.all
            print(repeat("\n", line_sep))
            print_remote_log(s, LOG_SERVER_INVERSE, input.log_lines)
            line_sep += skip_n_line - line_sep
        end
        if input.log_result || input.all
            print(repeat("\n", line_sep))
            print_remote_log(s, LOG_SERVER_RESULT, input.log_lines)
            line_sep += skip_n_line - line_sep
        end
    end
end

function print_log(logfile::String, nlines::Integer=5, indent::Integer=4)
    if !isfile(logfile)
        @warn "File '$logfile' not found"
        return nothing
    end
    buffer = readlines(logfile)
    L = length(buffer)
    indentstr = repeat(" ", indent)
    println(logfile, ":")
    foreach(buffer[max(1, L-nlines+1):end]) do l
        if contains(l, "INFO")
            printstyled(indentstr, l, "\n"; color=:light_green)
        elseif contains(l, "WARN")
            printstyled(indentstr, l, "\n"; color=:light_yellow)
        elseif contains(l, "ERROR")
            printstyled(indentstr, l, "\n"; color=:light_red)
        else
            println(indentstr, l)
        end
    end
end

function print_remote_log(s::InvServer, logfile::Function, nlines::Integer=5, indent::Integer=4)
    cmd = `ssh $(s.user)@$(s.ip) "cat $(logfile(s))"`

    texts = ""
    try
        texts = readchomp(cmd)
    catch err
        @error err
        return nothing
    end
    buffer = split(texts, '\n')
    indentstr = repeat(" ", indent)
    L = length(buffer)
    println(s.hostname, " ", logfile(s))
    foreach(buffer[max(1, L-nlines+1):end]) do l
        if contains(l, "INFO")
            printstyled(indentstr, l, "\n"; color=:light_green)
        elseif contains(l, "WARN")
            printstyled(indentstr, l, "\n"; color=:light_yellow)
        elseif contains(l, "ERROR")
            printstyled(indentstr, l, "\n"; color=:light_red)
        else
            println(indentstr, l)
        end
    end
end

function print_table(
    t::AbstractMatrix{<:AbstractString};
    delimiter::AbstractString=" ",
    head::Bool=false,
    indent::Integer=4,
    align::Union{Symbol,Vector{Symbol}}=:left,
    alignleft::AbstractVector{<:Integer} = Int[],
    aligncenter::AbstractVector{<:Integer} = Int[],
    alignright::AbstractVector{<:Integer} = Int[])
    if size(t, 1) == 0 || size(t, 2) == 0
        return nothing
    end
    column_width = vec(maximum(map(length, t), dims=1))
    if typeof(align) == Symbol
        _align = fill(align, size(column_width))
    else
        _align = align
    end
    _align[alignleft] .= :left
    _align[aligncenter] .= :center
    _align[alignright] .= :right
    buffer = map(CartesianIndices(t)) do I
        prefix = ""
        suffix = ""
        L = column_width[I[2]]-length(t[I])
        if _align[I[2]] == :left
            suffix = repeat(" ", L)
        elseif _align[I[2]] == :right
            prefix = repeat(" ", L)
        elseif _align[I[2]] == :center
            prefix = repeat(" ", floor(Int, L/2))
            suffix = repeat(" ", L - length(prefix))
        else
            suffix = repeat(" ", L)
        end
        return join([prefix, t[I], suffix])
    end
    indentstr = repeat(" ", indent)
    blines = map(l->indentstr*join(l, delimiter), eachrow(buffer))
    println(blines[1])
    if length(blines) == 1
        return nothing
    end
    if head
        println(indentstr, repeat("-", length(blines[1])))
    end
    foreach(println, blines[2:end])
    return nothing
end

try
    get_single_process_lock(@__FILE__)
    main()
catch
finally
    release_single_process_lock(@__FILE__)
end