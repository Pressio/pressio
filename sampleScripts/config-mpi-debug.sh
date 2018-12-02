#!/bin/bash

EXTRA_ARGS=$@
SRC=/Users/fnrizzi/Desktop/work/ROM/codes/sources/rompp
PFX=/Users/fnrizzi/Desktop/work/ROM/codes/installs/rompp_install

MPIPATH=/Users/fnrizzi/tpl/openmpi/301/installgcc730
TRILPATH=/Users/fnrizzi/tpl/trilinos/install_dyn_debug_gcc730_ompi301_static
EIGENINCPATH=/Users/fnrizzi/tpl/eigen/3.3.5/install
GTESTPATH=/Users/fnrizzi/tpl/gtest/installgcc730
BLAZEINCPATH=/Users/fnrizzi/tpl/blaze/3.4/install/include
ARMADILLOPATH=/Users/fnrizzi/tpl/armadillo/install_gcc730

cmake \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_INSTALL_PREFIX:PATH=${PFX} \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D TPL_FIND_SHARED_LIBS=OFF \
    -D rompp_LINK_SEARCH_START_STATIC=ON \
    \
    -D BLAS_LIBRARY_NAMES:STRING="openblas" \
    -D BLAS_LIBRARY_DIRS:PATH=/opt/local/lib \
    -D LAPACK_LIBRARY_NAMES:STRING="openblas" \
    -D LAPACK_LIBRARY_DIRS:PATH=/opt/local/lib \
    \
    -D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
    -D rompp_ENABLE_CXX11:BOOL=ON \
    -D rompp_ENABLE_SHADOW_WARNINGS:BOOL=OFF \
    -D CMAKE_CXX_FLAGS="-fopenmp" \
    \
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
    -D TPL_ENABLE_ARMADILLO=OFF \
    -D ARMADILLO_INCLUDE_DIRS:PATH=${ARMADILLOPATH}/include \
    -D ARMADILLO_LIBRARY_DIRS:PATH=${ARMADILLOPATH}/lib \
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
    -D rompp_ENABLE_qr:BOOL=ON \
    -D rompp_ENABLE_solvers:BOOL=ON \
    -D rompp_ENABLE_svd:BOOL=ON \
    -D rompp_ENABLE_ode:BOOL=ON \
    -D rompp_ENABLE_rom:BOOL=ON \
    \
    $EXTRA_ARGS \
    ${SRC}
