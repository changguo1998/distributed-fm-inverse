#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export JULIA_HOME="$DIR/julia"
export PATH="$JULIA_HOME/bin:$PATH"
export JULIA_DEPOT_PATH="$DIR/julia"
export JULIA_NUM_PRECOMPILE_TASKS=16

# echo "julia bin path: $(which julia)"
# julia --startup-file=no -E '
# function checkexist(K)
#     print("check $K exists: ")
#     if haskey(ENV, K)
#         printstyled(true; color=:light_green)
#         println("\t", ENV[K])
#     else
#         printstyled(false, "\n"; color=:red)
#     end
#     return haskey(ENV, K)
# end

# checkexist("JULIA_PKG_SERVER")

# if checkexist("JULIA_HOME") &&
#     checkexist("JULIA_DEPOT_PATH")
#     exit(0)
# else
#     exit(1)
# end
# '

julia_script="$1"

shift

if [ -f "$julia_script" ]; then
    exec julia --startup-file=no "$julia_script" "$@"
fi