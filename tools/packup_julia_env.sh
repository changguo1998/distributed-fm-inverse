#!/usr/bin/bash
set -eu
julia_bin_url="https://mirrors.ustc.edu.cn/julia-releases/bin/linux/x64/1.12/julia-1.12.2-linux-x86_64.tar.gz"
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

echo "root: $root"

if [ -d $root/julia ]; then
    rm -rf $root/julia
    mkdir -p $root/julia
fi

cd $root/tools/
if [ ! -f "julia.tar.gz" ]; then
    echo "Download julia interpreter"
    wget $julia_bin_url -O julia.tar.gz
fi

echo "Unpack julia interpreter"
mkdir -p tmp/
tar xaf julia.tar.gz -C tmp/
julia_bin_name="$(ls tmp/ | awk '{if(NR==1) print $0}')"
for f in `ls tmp/$julia_bin_name`; do
    mv tmp/$julia_bin_name/$f $root/julia
done
rm -rf tmp/

echo "Check environment"
cd $root
export JULIA_HOME=$root/julia
export PATH=$JULIA_HOME/bin:$PATH
export JULIA_DEPOT_PATH=$root/julia
export JULIA_PKG_SERVER="https://mirrors.cernet.edu.cn/julia"
echo "Bash env: JULIA_HOME: $JULIA_HOME"
echo "Bash env: JULIA_DEPOT_PATH: $JULIA_DEPOT_PATH"
echo "Bash env: JULIA_PKG_SERVER: $JULIA_PKG_SERVER"
echo "julia bin path: $(which julia)"
julia --startup-file=no -E '
function checkexist(K)
    print("check $K exists: ")
    if haskey(ENV, K)
        printstyled(true; color=:light_green)
        println("\t", ENV[K])
    else
        printstyled(false, "\n"; color=:red)
    end
    return haskey(ENV, K)
end

if checkexist("JULIA_HOME") &&
    checkexist("JULIA_DEPOT_PATH") &&
    checkexist("JULIA_PKG_SERVER")
    exit(0)
else
    exit(1)
end
'

echo "Install packages"
julia -t 32 --startup-file=no -e 'using Pkg; Pkg.activate("."; io=devnull); Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()'
if [ -f "env.tar.gz" ]; then
    rm env.tar.gz
fi

echo "Packaging Julia environment"
tar czf env.tar.gz julia
