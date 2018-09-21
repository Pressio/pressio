
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
#To include "rom_config.h" (because it lives in the binary folder)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/rom)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB rom_gen_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/rom/*.hpp")
FILE(GLOB rom_exp_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/rom/experimental/*.hpp")

FILE(GLOB rom_gen_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/rom/src/*.cc")
FILE(GLOB rom_exp_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/rom/experimental/src/*.cc")
    
SET(SOURCES 
  ${rom_gen_SRCS}
  ${rom_exp_SRCS}
)


SET(ROM_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

# #--------------------------------------
# SET(LUMPED_HEADERS
#   CORE_ALL
#   CORE_BASIC
#   CORE_MATRIX
#   CORE_MULTI_VECTOR
#   CORE_OPS
#   CORE_VECTOR
# )
# INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${ROM_INCDIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ROM_INCDIR}/rom)

INSTALL(FILES ${rom_gen_HEADERS} DESTINATION ${ROM_INCDIR}/rom)
INSTALL(FILES ${rom_exp_HEADERS} DESTINATION ${ROM_INCDIR}/rom/experimental)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  rom
  SOURCES ${SOURCES}
  )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
