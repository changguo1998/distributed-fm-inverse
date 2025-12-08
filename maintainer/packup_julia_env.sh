#!/usr/bin/bash
set -eu
julia_bin_url="https://mirrors.ustc.edu.cn/julia-releases/bin/linux/x64/1.12/julia-1.12.2-linux-x86_64.tar.gz"
cd ..
root="$(pwd)"
julia_env_path="$root/julia"

cd $root/maintainer/
mkdir -p tmp/
if [ -d $julia_env_path ]; then
    rm -rf $julia_env_path
    mkdir -p $julia_env_path
fi

if [ ! -f "julia.tar.gz" ]; then
    wget $julia_bin_url -O julia.tar.gz
fi
tar xaf julia.tar.gz -C tmp/
julia_bin_name="$(ls tmp/ | awk '{if(NR==1) print $0}')"
for f in `ls tmp/$julia_bin_name`; do
    mv tmp/$julia_bin_name/$f $julia_env_path
done
rm -rf tmp/

cd $root
JULIA_DEPOT_PATH=$julia_env_path
JULIA_PKG_SERVER="https://mirrors.cernet.edu.cn/julia"
$julia_env_path/bin/julia -t 8 --startup-file=no -e 'using Pkg; Pkg.activate("."); Pkg.instantiate();'

tar czf env.tar.gz julia
