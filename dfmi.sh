#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PATH="$DIR/julia/bin:$PATH"
JULIA_DEPOT_PATH="$DIR/julia"

julia_script="$1"

shift

if [ -f "$julia_script" ]; then
    exec julia --startup-file=no "$julia_script" "$@"
fi