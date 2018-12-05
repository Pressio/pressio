#!/bin/bash

EXTRA_ARGS=$@

SRC=/projects/rompp/sources/rompp
PFX=/projects/rompp/installs/rompp_install

MPIPATH=/projects/sierra/linux_rh6/SDK/mpi/openmpi/1.10.2-gcc-7.2.0-RHEL6/bin
TRILPATH=/gpfs1/fpierce/Trilinos-repo/install
EIGENINCPATH=/projects/rompp/tpls/eigen/3.3.5/install
GTESTPATH=/projects/rompp/tpls/gtest/install-sierra

cmake \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_INSTALL_PREFIX:PATH=${PFX} \
    \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D TPL_FIND_SHARED_LIBS=OFF \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
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
    -D TPL_ENABLE_GTEST=ON \
    -D GTEST_LIBRARY_DIRS:PATH=${GTESTPATH}/lib64 \
    -D GTEST_INCLUDE_DIRS:PATH=${GTESTPATH}/include \
    \
    -D rompp_ENABLE_Fortran=OFF \
    -D rompp_ENABLE_TESTS:BOOL=ON \
    -D rompp_ENABLE_EXAMPLES:BOOL=OFF \
    -D rompp_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D rompp_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    \
    -D rompp_ENABLE_core:BOOL=ON \
    -D rompp_ENABLE_qr:BOOL=OFF \
    -D rompp_ENABLE_solvers:BOOL=OFF \
    -D rompp_ENABLE_svd:BOOL=OFF \
    -D rompp_ENABLE_ode:BOOL=OFF \
    -D rompp_ENABLE_rom:BOOL=OFF \
    \
    $EXTRA_ARGS \
    ${SRC}