
option(TPL_ENABLE_ARMADILLO "Enable Armadillo TPL" OFF)
if(TPL_ENABLE_ARMADILLO)
  set(HAVE_ARMADILLO ON)
  set(ARMADILLO_LIB_NAMES armadillo)  
endif()
