#!/usr/bin/env bash

# make sure you have guardonce installed: https://github.com/cgmb/guardonce

# this script must be run where it is
# no args needed

# rel dir
PCK_REL_DIR=${PWD}/../include/pressio

# array of names
declare -a pcks=("mpl" "utils" "type_traits" "expressions")
# "ops" "qr" "solvers_linear" "solvers_nonlinear" "ode_advancers" "ode_steppers" "rom")
# loop over and fix he
for packName in ${pcks[@]}; do
    # target dir
    PCK_DIR=${PCK_REL_DIR}/${packName}
    echo ${PCK_DIR}

    cd ${PCK_DIR}
    echo ${PCK_DIR}

    # first, convert all guards to pragmas
    guard2once -r .
    # then converts from pragmas to header with specific pattern
    once2guard -r -p "path -2 | prepend ${packName}_ | upper | append _" -s "#endif  // %\n" .
done
