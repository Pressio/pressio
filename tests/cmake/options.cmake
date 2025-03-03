option(PRESSIO_ENABLE_CXX20 "Enable C++20" OFF)

# c++ standard
#=====================================================================
if (NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

message("${Magenta}>> PRESSIO: CMAKE_CXX_STANDARD is set to ${CMAKE_CXX_STANDARD} ${ColourReset}")

if(CMAKE_CXX_STANDARD STREQUAL "17"
    OR CMAKE_CXX_STANDARD STREQUAL "20")
  if(NOT "cxx_std_${CMAKE_CXX_STANDARD}" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(FATAL_ERROR "Compiler does not support C++${CMAKE_CXX_STANDARD}. This is required.")
  endif()
endif()

if(CMAKE_CXX_STANDARD EQUAL 20)
  add_definitions(-DPRESSIO_ENABLE_CXX20)
endif()

# pull in pressio-ops and read macros
#=====================================================================
if(PRESSIO_OPS_INCLUDE_DIR)
  message("${Green}-- Using ${PRESSIO_OPS_INCLUDE_DIR}")
else()
  include(FetchContent)
  message("${Green}-- PRESSIO_OPS_INCLUDE_DIR is not defined, so fetching pressio-ops...")
  FetchContent_Declare(
    pressio-ops
    GIT_REPOSITORY https://github.com/Pressio/pressio-ops.git
    GIT_TAG        develop
    GIT_PROGRESS   TRUE
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(pressio-ops)
  set(PRESSIO_OPS_INCLUDE_DIR ${pressio-ops_SOURCE_DIR}/include)
endif()
include_directories(SYSTEM ${PRESSIO_OPS_INCLUDE_DIR})

# logging
#=====================================================================
option(PRESSIO_ENABLE_LOGGING "Enable logging via pressio-log" ON)
if (PRESSIO_ENABLE_LOGGING)
  add_compile_definitions(PRESSIO_ENABLE_LOGGING=1)
else()
  add_compile_definitions(PRESSIO_ENABLE_LOGGING=0)
endif()

# get pressio-log (either find it or fetch it)
message(STATUS "PRESSIO_ENABLE_LOGGING = ${PRESSIO_ENABLE_LOGGING}")
if (PRESSIO_ENABLE_LOGGING)
  if(PRESSIO_LOG_INCLUDE_DIR)
    message("${Green}-- Using ${PRESSIO_LOG_INCLUDE_DIR}")
  else()
    include(FetchContent)
    message("${Green}-- PRESSIO_LOG_INCLUDE_DIR is not defined, so fetching pressio-log...")
    FetchContent_Declare(
      pressio-log
      GIT_REPOSITORY https://github.com/Pressio/pressio-log.git
      GIT_TAG        main
      GIT_PROGRESS   TRUE
      GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(pressio-log)
    set(PRESSIO_LOG_INCLUDE_DIR ${pressio-log_SOURCE_DIR}/include)
  endif()
  include_directories(SYSTEM ${PRESSIO_LOG_INCLUDE_DIR})
else()
  message(STATUS "Logging has been disabled (PRESSIO_ENABLE_LOGGING=OFF)")
endif()

# TPLs
#=====================================================================
option(PRESSIO_ENABLE_TPL_EIGEN		  "Enable Eigen TPL"	  ON)
option(PRESSIO_ENABLE_TPL_TRILINOS	"Enable Trilinos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_KOKKOS		"Enable Kokkos TPL"	  OFF)
option(PRESSIO_ENABLE_TPL_MPI		    "Enable MPI"	      	OFF)

if(PRESSIO_ENABLE_TPL_EIGEN)
  message(">> Eigen is currently enabled by default via PRESSIO_ENABLE_TPL_EIGEN=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_EIGEN)

  if(EIGEN_INCLUDE_DIR)
    include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})
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

  include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
  link_libraries(${Trilinos_LIBRARIES})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  message(">> Enabling Kokkos since PRESSIO_ENABLE_TPL_KOKKOS=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_KOKKOS)

  # when trilinos is enabled it links kokkos too
  if(NOT PRESSIO_ENABLE_TPL_TRILINOS)
    # if kokkos is used as standalone lib, then we are more specific
    # user needs to define: Kokkos_ROOT and KokkosKernels_ROOT
    if (NOT Kokkos_ROOT OR NOT KokkosKernels_ROOT)
      message(
        FATAL_ERROR
        "Missing Kokkos_ROOT. Kokkos needs:
        -D Kokkos_ROOT=<full-path-to-kokkos-installation>
        -D KokkosKernels_ROOT=<full-path-to-kokkos-kernels-installation>"
      )
    endif()
    find_package(Kokkos REQUIRED)
    find_package(KokkosKernels REQUIRED)
    link_libraries(Kokkos::kokkos Kokkos::kokkoskernels)
  endif()
endif()

if(PRESSIO_ENABLE_TPL_MPI)
  message(">> Enabling MPI since PRESSIO_ENABLE_TPL_MPI=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_MPI)

  find_package(MPI REQUIRED)
  include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
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
