
option(TPL_ENABLE_TRILINOS "Enable Trilinos TPL" OFF)

if(TPL_ENABLE_TRILINOS)
  # this define is used by the cmake_config
  set(HAVE_TRILINOS ON)

  # if we need to build tests, then add include and
  if(BUILD_UNIT_TESTS OR BUILD_TESTS)

    set(HAVE_TEUCHOS_TIMERS ON)
    if( (NOT TRILINOS_INC_DIR OR NOT TRILINOS_LIB_DIR)
	AND (NOT TRILINOS_INCLUDE_DIR OR NOT TRILINOS_LIBRARIES_DIR))
      message(FATAL_ERROR
	"TPL_ENABLE_TRILINOS=${TPL_ENABLE_TRILINOS} but did not specify how to find it.
        Please reconfigure with:
          -DTRILINOS_INC_DIR=<full-path-to-headers>
          -DTRILINOS_LIB_DIR=<full-path-to-libs>
          or
          -TRILINOS_INCLUDE_DIR=<full-path-to-headers>
          -TRILINOS_LIBRARIES_DIR=<full-path-to-libs>
          ")
    endif()

    if(NOT TRILINOS_INC_DIR AND TRILINOS_INCLUDE_DIR)
      set(TRILINOS_INC_DIR ${TRILINOS_INCLUDE_DIR})
    endif()

    if(NOT TRILINOS_LIB_DIR AND TRILINOS_LIBRARIES_DIR)
      set(TRILINOS_LIB_DIR ${TRILINOS_LIBRARIES_DIR})
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

    include_directories(${TRILINOS_INC_DIR})
    link_directories(${TRILINOS_LIB_DIR})
    link_libraries(${TRILINOS_LIB_NAMES})
  endif()
endif()
