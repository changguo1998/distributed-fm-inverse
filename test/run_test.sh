#!/usr/bin/bash
root=/home/guochang/Projects/distributed-fm-inverse/
set -e
N=4

reset
clear
echo "Reset test environment"
cd $root/test/
julia cleanup_test_temporary_file.jl
julia generate_test_event.jl $N

cd $root
for((i=0; i<=$N; i++)); do
    bash dfmi.sh host/submit_to_host_buffer.jl
done
bash dfmi.sh server/update_server_status.jl
for((i=0; i<=$N; i++)); do
    bash dfmi.sh host/upload_to_server.jl
done
bash dfmi.sh server/update_server_status.jl

exit 0
set -x

bash dfmi.sh server/unpack_input_file.jl
bash dfmi.sh server/call_inverse.jl
bash dfmi.sh server/pack_result.jl
