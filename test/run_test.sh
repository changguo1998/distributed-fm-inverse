set -x
N=3

echo "Reset test environment"
cd /home/guochang/Projects/distributed-fm-inverse/
julia reset.jl

echo "Test server"
cd /home/guochang/Projects/distributed-fm-inverse/server
julia update_server_status.jl
julia unpack_input_file.jl
julia call_inverse.jl
julia pack_result.jl
exit(0)

echo "Test host"
cd /home/guochang/Projects/distributed-fm-inverse/host

for((i=0; i<$N; i++)); do
    julia submit_to_host_buffer.jl
done

julia upload_to_server.jl
