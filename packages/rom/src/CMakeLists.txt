
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "core_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

FILE(GLOB rom_gen_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB rom_exp_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.hpp")

FILE(GLOB rom_gen_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
FILE(GLOB rom_exp_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/src/*.cc")
    
SET(SOURCES 
  ${rom_gen_SRCS}
  ${rom_exp_SRCS}
)


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
SET(ROM_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/rom/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ROM_INCDIR})
INSTALL(FILES ${rom_gen_HEADERS} DESTINATION ${ROM_INCDIR})
INSTALL(FILES ${rom_exp_HEADERS} DESTINATION ${ROM_INCDIR}/experimental)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  rom
  #   HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
