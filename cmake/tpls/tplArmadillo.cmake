
option(PRESSIO_ENABLE_TPL_ARMADILLO "Enable Armadillo TPL" OFF)

if(PRESSIO_ENABLE_TPL_ARMADILLO)

  # if we need to build tests, then prep for it
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    set(ARMADILLO_LIB_NAMES armadillo)
  endif()
endif()
