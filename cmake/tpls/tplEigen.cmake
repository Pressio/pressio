

if(PRESSIO_ENABLE_TPL_EIGEN)
  message(">> Eigen is currently enabled by default via PRESSIO_ENABLE_TPL_EIGEN=ON")
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
