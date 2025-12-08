#!/usr/bin/bash
root=/home/guochang/Projects/distributed-fm-inverse/
set -e
set -x
N=3

echo "Reset test environment"
cd $root
julia reset.jl

cd $root/test/
julia generate_test_event.jl $N

cd $root/host
for((i=0; i<$N; i++)); do
    julia submit_to_host_buffer.jl
done

cd $root/server
julia update_server_status.jl

cd $root/host
julia upload_to_server.jl

cd $root/server
julia update_server_status.jl
julia unpack_input_file.jl
julia call_inverse.jl
julia pack_result.jl
exit 0
