
option(PRESSIO_ENABLE_TPL_PYBIND11 "Enable Pybind11 TPL" OFF)
if(PRESSIO_ENABLE_TPL_PYBIND11)
  message("Found PRESSIO_ENABLE_TPL_PYBIND11=${PRESSIO_ENABLE_TPL_PYBIND11}.")

  # if we need to build tests, then prep for it
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    if(NOT PYBIND11_DIR AND NOT PYBIND11_INCLUDE_DIR)
      message(FATAL_ERROR
      "I cannot find the Pybind11 headers. Please reconfigure with:
      	-DPYBIND11_DIR=<full-path-to-toplevel-install>
      	-DPYBIND11_INCLUDE_DIR=<full-path-to-headers>
      ")
    endif()

    find_package(pybind11 REQUIRED PATHS ${PYBIND11_DIR}/share/cmake)
    find_package(Python3 COMPONENTS Interpreter NumPy)
    if(NOT PYBIND11_INC_DIR AND PYBIND11_INCLUDE_DIR)
      set(PYBIND11_INC_DIR ${PYBIND11_INCLUDE_DIR})
    endif()
    include_directories(${Python3_INCLUDE_DIRS} ${PYBIND11_INC_DIR})

    if(NOT ${Python3_FOUND})
      message(FATAL_ERROR "Python > 3 not found")
    endif()
  endif()
endif()
