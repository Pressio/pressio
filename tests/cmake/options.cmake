option(PRESSIO_ENABLE_CXX14 "Enable C++14" OFF)
option(PRESSIO_ENABLE_CXX17 "Enable C++17" OFF)
option(PRESSIO_ENABLE_CXX20 "Enable C++20" OFF)

if (CMAKE_CXX_STANDARD EQUAL 14)
  add_definitions(-DPRESSIO_ENABLE_CXX14)
elseif(CMAKE_CXX_STANDARD EQUAL 17)
  add_definitions(-DPRESSIO_ENABLE_CXX17)
elseif(CMAKE_CXX_STANDARD EQUAL 20)
  add_definitions(-DPRESSIO_ENABLE_CXX20)
endif()

option(PRESSIO_ENABLE_DEBUG_PRINT "Enable debug printing" OFF)
if (PRESSIO_ENABLE_DEBUG_PRINT)
  add_definitions(-DPRESSIO_ENABLE_DEBUG_PRINT)
endif()

option(PRESSIO_ENABLE_TPL_EIGEN		  "Enable Eigen TPL"	  ON)
option(PRESSIO_ENABLE_TPL_TRILINOS	"Enable Trilinos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_KOKKOS		"Enable Kokkos TPL"	  OFF)
option(PRESSIO_ENABLE_TPL_MPI		    "Enable MPI"	      	OFF)
option(PRESSIO_ENABLE_TPL_PYBIND11	"Enable Pybind11 TPL"	OFF)


if(PRESSIO_ENABLE_TPL_EIGEN)
  message(">> Eigen is currently enabled by default via PRESSIO_ENABLE_TPL_EIGEN=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_EIGEN)

  if(EIGEN_INCLUDE_DIR)
    include_directories(${EIGEN_INCLUDE_DIR})
  else()
    find_package(Eigen3)

    if(NOT EIGEN3_FOUND)
      # TODO: use FetchContent_Declare instead of failing?
      message(FATAL_ERROR
      "I cannot find the Eigen headers. "
      "Please reconfigure with: -DEIGEN_INCLUDE_DIR=<full-path-to-headers>")
    endif()

    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
  endif()
endif()

# if trilinos is on, then also set MPI, BLAS, LAPACK and KOKKOS ON
if(PRESSIO_ENABLE_TPL_TRILINOS)
  message(">> PRESSIO_ENABLE_TPL_TRILINOS=ON ==> enabling also BLAS, LAPACK, MPI, KOKKOS")
  add_definitions(-DPRESSIO_ENABLE_TPL_TRILINOS)

  set(PRESSIO_ENABLE_TPL_KOKKOS ON)
  set(PRESSIO_ENABLE_TPL_MPI ON)
  set(PRESSIO_ENABLE_TPL_BLAS ON)
  set(PRESSIO_ENABLE_TPL_LAPACK ON)

  find_package(Trilinos REQUIRED)
  # TODO: it is possible to use find_package(<PackageName>) for each (sub)package
  # https://trilinos.github.io/pdfs/Finding_Trilinos.txt

  include_directories(${Trilinos_INCLUDE_DIRS})
  link_libraries(${Trilinos_LIBRARIES})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  message(">> Enabling Kokkos since PRESSIO_ENABLE_TPL_KOKKOS=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_KOKKOS)

  # when trilinos is enabled it links kokkos too
  if(NOT PRESSIO_ENABLE_TPL_TRILINOS)
    # if kokkos is used as standalone lib, then we are more specific
    # user needs to define: KOKKOS_ROOT_DIR and KOKKOS_KERNELS_ROOT_DIR
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

if(PRESSIO_ENABLE_TPL_MPI)
  message(">> Enabling MPI since PRESSIO_ENABLE_TPL_MPI=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_MPI)

  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_DIRS})
endif()

if(PRESSIO_ENABLE_TPL_Pybind11)
  message(">> Enabling Pybind11 since PRESSIO_ENABLE_TPL_PYBIND11=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_Pybind11)
endif()

if(PRESSIO_ENABLE_TPL_BLAS)
  message(">> Enabling BLAS since PRESSIO_ENABLE_TPL_BLAS=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_BLAS)

  cmake_policy(SET CMP0074 NEW)
  find_package(BLAS REQUIRED)
  link_libraries(${BLAS_LIBRARIES})
  message("BLASLIBS=${BLAS_LIBRARIES}")
endif()

if(PRESSIO_ENABLE_TPL_LAPACK)
  message(">> Enabling LAPACK since PRESSIO_ENABLE_TPL_LAPACK=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_LAPACK)

  cmake_policy(SET CMP0074 NEW)
  find_package(LAPACK REQUIRED)
  link_libraries(${LAPACK_LIBRARIES})
  message("LAPLIBS=${LAPACK_LIBRARIES}")
endif()
