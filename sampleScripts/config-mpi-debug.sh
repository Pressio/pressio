#!/bin/bash

EXTRA_ARGS=$@

SRC=/Users/fnrizzi/Desktop/cleanup/codes/rompp
PFX=/Users/fnrizzi/Desktop/romppInstall

MPIPATH=/Users/fnrizzi/tpl/openmpi/301/installgcc640
TRILPATH=/Users/fnrizzi/tpl/trilinos/install_debug_gcc640_ompi301
EIGENINCPATH=/Users/fnrizzi/tpl/eigen/3.3.5/install
BLAZEINCPATH=/Users/fnrizzi/tpl/blaze/3.4/install/include
ARMADILLOPATH=/Users/fnrizzi/tpl/armadillo/install_gcc640
KOKKOSPATH=/Users/fnrizzi/tpl/kokkos/install
GTESTPATH=/Users/fnrizzi/tpl/gtest/install

cmake \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_INSTALL_PREFIX:PATH=${PFX} \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    \
    -D MPI_EXEC_MAX_NUMPROCS:STRING=16 \
    -D rompp_ENABLE_CXX11:BOOL=ON\
    -D rompp_ENABLE_SHADOW_WARNINGS:BOOL=OFF \
    \
    -D CMAKE_CXX_FLAGS="-std=c++14 -fopenmp" \
    -D TPL_ENABLE_MPI=ON \
    -D MPI_BASE_DIR:PATH=${MPIPATH} \
    -D MPI_EXEC:FILEPATH=${MPIPATH}/bin/mpirun \
    -D TPL_ENABLE_TRILINOS=ON \
    -D TRILINOS_LIBRARY_DIRS:PATH=${TRILPATH}/lib \
    -D TRILINOS_INCLUDE_DIRS:PATH=${TRILPATH}/include \
    -D TPL_ENABLE_EIGEN=ON \
    -D EIGEN_INCLUDE_DIRS:PATH=${EIGENINCPATH} \
    -D TPL_ENABLE_BLAZE=OFF \
    -D BLAZE_INCLUDE_DIRS:PATH=${BLAZEINCPATH} \
    -D TPL_ENABLE_ARMADILLO=ON \
    -D ARMADILLO_INCLUDE_DIRS:PATH=${ARMADILLOPATH}/include \
    -D ARMADILLO_LIBRARY_DIRS:PATH=${ARMADILLOPATH}/lib \
    \
    -D TPL_ENABLE_GTEST=ON \
    -D GTEST_LIBRARY_DIRS:PATH=${GTESTPATH}/lib \
    -D GTEST_INCLUDE_DIRS:PATH=${GTESTPATH}/include \
    \
    -D rompp_ENABLE_Fortran=OFF \
    -D rompp_ENABLE_TESTS:BOOL=ON \
    -D rompp_ENABLE_EXAMPLES:BOOL=OFF \
    -D rompp_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D rompp_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    \
    -D rompp_ENABLE_core:BOOL=ON \
    -D rompp_ENABLE_solvers:BOOL=ON \
    -D rompp_ENABLE_svd:BOOL=ON \
    -D rompp_ENABLE_ode:BOOL=ON \
    -D rompp_ENABLE_rom:BOOL=ON \
    \
    $EXTRA_ARGS \
    ${SRC}

    #-D rompp_ENABLE_optimization:BOOL=ON \
    #-D rompp_ENABLE_apps:BOOL=ON \
    # -DTribitsExProj_TRIBITS_DIR=$TRIBITS_DIR \
    # -D rompp_ENABLE_svd:BOOL=OFF \
    # -D rompp_ENABLE_apps:BOOL=OFF \
    # -D rompp_ENABLE_rom:BOOL=OFF \
    # -D TPL_ENABLE_KOKKOS=ON \
    # -D KOKKOS_LIBRARY_DIRS:PATH=${KOKKOSPATH}/lib \
    # -D KOKKOS_INCLUDE_DIRS:PATH=${KOKKOSPATH}/include \
    #-D rompp_ENABLE_CXX11:BOOL=ON\
    #-D CMAKE_CXX_FLAGS="-std=c++11 -fopenmp"\
