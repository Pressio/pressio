

# if trilinos is on, then also set MPI, BLAS, LAPACK and KOKKOS ON
if(PRESSIO_ENABLE_TPL_TRILINOS)
  message(">> PRESSIO_ENABLE_TPL_TRILINOS=ON ==> enabling also BLAS, LAPACK, MPI, KOKKOS")

  set(PRESSIO_ENABLE_TPL_KOKKOS ON)
  set(PRESSIO_ENABLE_TPL_MPI ON)
  set(PRESSIO_ENABLE_TPL_BLAS ON)
  set(PRESSIO_ENABLE_TPL_LAPACK ON)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  #message("Enabling Trilinos since PRESSIO_ENABLE_TPL_TRILINOS=${PRESSIO_ENABLE_TPL_TRILINOS}")

  # if we need to build tests, then find trilinos
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
