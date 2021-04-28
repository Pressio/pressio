

if(PRESSIO_ENABLE_TPL_Pybind11)
  message(">> Enabling Pybind11 since PRESSIO_ENABLE_TPL_PYBIND11=ON")
endif()

function(pybind11_fatal)
  message(FATAL_ERROR "I cannot find the Pybind11 library.
To use it, please reconfigure with:
-D PYBIND11_ROOT=<full-path-to-installation>")
endfunction()


if(PRESSIO_ENABLE_TPL_PYBIND11)
  message("Found PRESSIO_ENABLE_TPL_PYBIND11=${PRESSIO_ENABLE_TPL_PYBIND11}.")

  # if we need to build tests, then prep for it
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    # if PYBIND11_ROOT not found
    if (NOT PYBIND11_ROOT)
      pybind11_fatal()
    endif()

    find_package(pybind11 REQUIRED PATHS ${PYBIND11_ROOT}/share/cmake)
    find_package(Python3 COMPONENTS Interpreter NumPy)

    include_directories(${Python3_INCLUDE_DIRS} ${PYBIND11_ROOT}/include)

    if(NOT ${Python3_FOUND})
      message(FATAL_ERROR "Python > 3 not found")
    endif()
  endif()
endif()
