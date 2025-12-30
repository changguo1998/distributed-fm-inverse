#!/usr/bin/env julia
include(joinpath(@__DIR__, "../lib.jl"))

function _rm(f::AbstractString)
    rm(f; recursive=true, force=true)
    return nothing
end

function main()
    rm_list = filter(readdir(BUFFER_SERVER_CLEAN())) do f
        s = stat(joinpath(BUFFER_SERVER_CLEAN(), f))
        return now() - unix2datetime(max(s.mtime, s.ctime)) > Second(2)
    end

    foreach(rm_list) do f
        tag = replace(f, ".flag"=>"")
        waitinglist = String[]
        fn = joinpath(BUFFER_SERVER_INPUT(), tag * "_input.tar.gz")
        if isfile(fn)
            push!(waitinglist, fn)
        end
        fn = joinpath(BUFFER_SERVER_RUN(), tag)
        if isdir(fn)
            push!(waitinglist, fn)
        end
        fn = joinpath(BUFFER_SERVER_RESULT(), tag * "_result.tar.gz")
        if isfile(fn)
            push!(waitinglist, fn)
        end
        if isempty(waitinglist)
            _rm(joinpath(BUFFER_SERVER_CLEAN(), f))
        else
            foreach(_rm, waitinglist)
        end
    end
end


try
    get_single_process_lock(@__FILE__)
    main()
catch err
    @error err
finally
    release_single_process_lock(@__FILE__)
end
