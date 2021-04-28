#!/usr/bin/env bash

# make sure you have guardonce installed: https://github.com/cgmb/guardonce

# this script must be run where it is
# no args needed

# rel dir
PCK_REL_DIR=${PWD}/../include

# array of names
declare -a pcks=("mpl" "utils" "containers" "ops" "apps" "qr" "svd" "optimizers" "solvers" "ode" "rom")
# loop over and fix he
for packName in ${pcks[@]}; do
    # target dir
    PCK_DIR=${PCK_REL_DIR}/${packName}
    echo ${PCK_DIR}

    cd ${PCK_DIR}

    # first, convert all guards to pragmas
    guard2once -r ./src
    # then converts from pragmas to header with specific pattern
    once2guard -r -p "path -2 | prepend ${packName}_ | upper | append _" -s "#endif  // %\n" ./src/
done
