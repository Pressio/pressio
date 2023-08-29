# if(PRESSIO_ENABLE_TPL_BLAS)
#   message(">> Enabling BLAS since PRESSIO_ENABLE_TPL_BLAS=ON ==> enabling also LAPACK")
#   set(PRESSIO_ENABLE_TPL_LAPACK ON)
# endif()
# if(PRESSIO_ENABLE_TPL_LAPACK)
#   message(">> Enabling LAPACK since PRESSIO_ENABLE_TPL_LAPACK=ON ==> enabling also BLAS")
#   set(PRESSIO_ENABLE_TPL_BLAS ON)
# endif()

if(PRESSIO_ENABLE_TPL_MPI)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
  endif()
endif()


if(PRESSIO_ENABLE_TPL_BLAS OR PRESSIO_ENABLE_TPL_TRILINOS)
    # check if BLAS_ROOT is specified
    if (NOT ${BLAS_ROOT})
      message("")
      message(FATAL_ERROR "BLAS_ROOT not speificed, terminating")
      message("Make sure you set the BLAS_ROOT env var")
    endif()

    cmake_policy(SET CMP0074 NEW)
    find_package(BLAS REQUIRED)
    link_libraries(${BLAS_LIBRARIES})
    message("")
endif()

if(PRESSIO_ENABLE_TPL_LAPACK)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    # check if LAPACK_ROOT is specified
    if (NOT ${LAPACK_ROOT})
      message("")
      message(FATAL_ERROR "LAPACK_ROOT not speificed, terminating")
      message("Make sure you set the LAPACK_ROOT env var")
    endif()

    cmake_policy(SET CMP0074 NEW)
    find_package( LAPACK REQUIRED )
    link_libraries(${LAPACK_LIBRARIES})
    message("LAPLIBS=${LAPACK_LIBRARIES}")
  endif()
endif()


if(PRESSIO_ENABLE_TPL_EIGEN)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    if(NOT EIGEN_INC_DIR AND NOT EIGEN_INCLUDE_DIR)
      message(FATAL_ERROR
	"I cannot find the Eigen headers. Please reconfigure with:
      	-DEIGEN_INC_DIR=<full-path-to-headers>
      	or
      	-DEIGEN_INCLUDE_DIR=<full-path-to-headers>
      ")
    endif()

    if(NOT EIGEN_INC_DIR AND EIGEN_INCLUDE_DIR)
      set(EIGEN_INC_DIR ${EIGEN_INCLUDE_DIR})
    endif()

    include_directories(${EIGEN_INC_DIR})
  endif()
endif()


if(PRESSIO_ENABLE_TPL_KOKKOS)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    # # when trilinos is also enabled it links kokkos too, see tplTrilinos.cmake
    if(NOT PRESSIO_ENABLE_TPL_TRILINOS)
      # if kokkos is used as standalone lib, then we are more specific
      # user needs to defined: KOKKOS_ROOT_DIR and KOKKOS_KERNELS_ROOT_DIR
      if (NOT KOKKOS_ROOT OR NOT KOKKOS_KERNELS_ROOT)
	message(FATAL_ERROR "Missing KOKKOS_ROOT. KOKKOS needs:
          -D KOKKOS_ROOT=<full-path-to-kokkos-installation>
          -D KOKKOS_KERNELS_ROOT=<full-path-to-kokkos-kernels-installation>
          ")
      endif()

      set(KOKKOS_LIB_NAMES kokkoscontainers kokkoscore kokkoskernels)

      include_directories(${KOKKOS_ROOT}/include ${KOKKOS_KERNELS_ROOT}/include)

      link_directories(${KOKKOS_ROOT}/lib ${KOKKOS_ROOT}/lib64
	${KOKKOS_KERNELS_ROOT}/lib ${KOKKOS_KERNELS_ROOT}/lib64)

      link_libraries(${KOKKOS_LIB_NAMES})
    endif()

  endif()
endif()


if(PRESSIO_ENABLE_TPL_TRILINOS)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    if (NOT TRILINOS_ROOT)
      message(FATAL_ERROR
	"You enabled PRESSIO_ENABLE_TPL_TRILINOS but did not set TRILINOS_ROOT.
        Please reconfigure with:
          -D TRILINOS_ROOT=<full-path-to-trilinos-install>
          ")
    endif()

    set(TRILINOS_LIB_NAMES kokkosalgorithms
      kokkoscontainers
      kokkoscore
      teuchoskokkoscomm
      teuchoskokkoscompat
      teuchosremainder
      teuchosnumerics
      teuchoscomm
      teuchosparameterlist
      teuchosparser
      teuchoscore
      epetra
      epetraext
      ifpack
      aztecoo
      tpetraext
      tpetrainout
      tpetra
      kokkostsqr
      tpetraclassiclinalg
      tpetraclassicnodeapi
      tpetraclassic
      kokkoskernels
      ifpack2
      triutils
      # repeat to solve issue we have on linux
      kokkosalgorithms
      teuchosparameterlist)

    include_directories(${TRILINOS_ROOT}/include)
    link_directories(${TRILINOS_ROOT}/lib ${TRILINOS_ROOT}/lib64)
    link_libraries(${TRILINOS_LIB_NAMES})
  endif()
endif()


# function(pybind11_fatal)
#   message(FATAL_ERROR "I cannot find the Pybind11 library.
# To use it, please reconfigure with:
# -D PYBIND11_ROOT=<full-path-to-installation>")
# endfunction()

# if(PRESSIO_ENABLE_TPL_PYBIND11)
#   message("Found PRESSIO_ENABLE_TPL_PYBIND11=${PRESSIO_ENABLE_TPL_PYBIND11}.")

#   # if we need to build tests, then prep for it
#   if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

#     # if PYBIND11_ROOT not found
#     if (NOT PYBIND11_ROOT)
#       pybind11_fatal()
#     endif()

#     find_package(pybind11 REQUIRED PATHS ${PYBIND11_ROOT}/share/cmake)
#     find_package(Python3 COMPONENTS Interpreter NumPy)

#     include_directories(${Python3_INCLUDE_DIRS} ${PYBIND11_ROOT}/include)

#     if(NOT ${Python3_FOUND})
#       message(FATAL_ERROR "Python > 3 not found")
#     endif()
#   endif()
# endif()