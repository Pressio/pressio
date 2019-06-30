#!/bin/bash

#export PYBPATH=/Users/fnrizzi/Desktop/work/ROM/clang/clang700_ompi400_dbg_shared/pybind11/install/include
export MPIPATH=/Users/fnrizzi/tpl/openmpi/4.0.0/install_clang700/bin
export CXX=${MPIPATH}/mpic++
${CXX} -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` main.cc -o example`python3-config --extension-suffix`
