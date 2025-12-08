#!/usr/bin/bash

if [ ! -d "julia" ]; then
    tar xaf env.tar.gz
fi
bash dfmi.sh reset.jl