#!/bin/bash

EXTRA_ARGS=$@

SRC=/Users/fnrizzi/Desktop/indwork/NGA/ROM/codes/rompp
PFX=/Users/fnrizzi/Desktop/romppInstall

MPIPATH=/Users/fnrizzi/tpl/openmpi/301/installgcc550
TRILPATH=/Users/fnrizzi/tpl/trilinos/installCPPonly

cmake \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_INSTALL_PREFIX:PATH=${PFX} \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    \
    -D TPL_ENABLE_MPI=ON \
    -D MPI_BASE_DIR:PATH=${MPIPATH} \
    -D MPI_EXEC:FILEPATH=${MPIPATH}/bin/mpirun \
    -D TPL_ENABLE_TRILINOS=ON \
    -D TRILINOS_LIBRARY_DIRS:PATH=${TRILPATH}/lib \
    -D TRILINOS_INCLUDE_DIRS:PATH=${TRILPATH}/include \
    \
    -D rompp_ENABLE_Fortran=OFF \
    -D rompp_ENABLE_TESTS:BOOL=ON \
    -D rompp_ENABLE_EXAMPLES:BOOL=ON \
    -D rompp_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D rompp_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    \
    -D rompp_ENABLE_commonBase:BOOL=ON \
    \
    $EXTRA_ARGS \
    ${SRC}

#    -DTribitsExProj_TRIBITS_DIR=$TRIBITS_DIR \
