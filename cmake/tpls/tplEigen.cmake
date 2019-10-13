

option(TPL_ENABLE_EIGEN "Enable Eigen TPL" ON)

if(TPL_ENABLE_EIGEN)
  set(HAVE_EIGEN ON)
  message("Found TPL_ENABLE_EIGEN=${TPL_ENABLE_EIGEN}, so HAVE_EIGEN=${HAVE_EIGEN}")

  # if we need to build tests, then add include and
  if(BUILD_UNIT_TESTS OR BUILD_TESTS)

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
