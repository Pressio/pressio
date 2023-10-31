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

option(PRESSIO_ENABLE_TPL_EIGEN		"Enable Eigen TPL"	ON)
option(PRESSIO_ENABLE_TPL_TRILINOS	"Enable Trilinos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_KOKKOS		"Enable Kokkos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_MPI		"Enable MPI"		OFF)
option(PRESSIO_ENABLE_TPL_PYBIND11	"Enable Pybind11 TPL"	OFF)
# option(PRESSIO_ENABLE_TEUCHOS_TIMERS "bla bla" OFF)


if(PRESSIO_ENABLE_TPL_EIGEN)
  message(">> Eigen is currently enabled by default via PRESSIO_ENABLE_TPL_EIGEN=ON")
  add_definitions(-DPRESSIO_ENABLE_TPL_EIGEN)

  if(EIGEN_INCLUDE_DIR)
    include_directories(${EIGEN_INCLUDE_DIR})
  else()
    find_package(Eigen3)

    if(NOT EIGEN3_FOUND)
      # FIXME: use FetchContent_Declare instead of failing?
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

  set(PRESSIO_ENABLE_TPL_KOKKOS ON)
  set(PRESSIO_ENABLE_TPL_MPI ON)

  find_package(Trilinos REQUIRED)
  # TODO: it is possible to use find_package(<PackageName>) for each (sub)package
  # https://trilinos.github.io/pdfs/Finding_Trilinos.txt
  add_definitions(-DPRESSIO_ENABLE_TPL_TRILINOS)

  include_directories(${Trilinos_INCLUDE_DIRS})
  link_libraries(${Trilinos_LIBRARIES})

  # FINISH THIS
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  message(">> Enabling Kokkos since PRESSIO_ENABLE_TPL_KOKKOS=ON")
endif()

if(PRESSIO_ENABLE_TPL_MPI)
  message(">> Enabling MPI since PRESSIO_ENABLE_TPL_MPI=ON")

  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_DIRS})
endif()

if(PRESSIO_ENABLE_TPL_Pybind11)
  message(">> Enabling Pybind11 since PRESSIO_ENABLE_TPL_PYBIND11=ON")
endif()


# if(PRESSIO_ENABLE_TPL_KOKKOS AND NOT PRESSIO_ENABLE_TPL_TRILINOS)
#   find_package(KokkosKernels REQUIRED)
#   set(KOKKOS_LIBS Kokkos::kokkoskernels)
#   #Kokkos::BLAS Kokkos::LAPACK Kokkos::kokkos Kokkos::kokkoskernels)
# endif()
