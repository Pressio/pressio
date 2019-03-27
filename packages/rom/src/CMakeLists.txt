
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "core_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

FILE(GLOB rom_gen_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

FILE(GLOB rom_decoder_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/decoder/*.hpp")
FILE(GLOB rom_decorators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/decorators/*.hpp")
FILE(GLOB rom_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/policies/*.hpp")

FILE(GLOB rom_galerkin_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/galerkin/*.hpp")
FILE(GLOB rom_lspg_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/*.hpp")
FILE(GLOB rom_lspg_generator_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/generator/*.hpp")

FILE(GLOB rom_lspg_steady_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/*.hpp")
FILE(GLOB rom_lspg_steady_generator_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/generator/*.hpp")

# sources
FILE(GLOB rom_gen_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
SET(SOURCES ${rom_gen_SRCS})

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ROM_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ROM_GALERKIN
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ROM_LSPG
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ROM_LSPG_STEADY
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
SET(ROM_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/rom/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ROM_INCDIR})
INSTALL(FILES ${rom_gen_HEADERS} DESTINATION ${ROM_INCDIR})

INSTALL(FILES ${rom_decoder_HEADERS} DESTINATION ${ROM_INCDIR}/decoder)
INSTALL(FILES ${rom_decorators_HEADERS} DESTINATION ${ROM_INCDIR}/decorators)
INSTALL(FILES ${rom_policies_HEADERS} DESTINATION ${ROM_INCDIR}/policies)

INSTALL(FILES ${rom_galerkin_HEADERS} DESTINATION ${ROM_INCDIR}/galerkin)

INSTALL(FILES ${rom_lspg_HEADERS} DESTINATION ${ROM_INCDIR}/lspg)
INSTALL(FILES ${rom_lspg_generator_HEADERS} DESTINATION ${ROM_INCDIR}/lspg/generator)

INSTALL(FILES ${rom_lspg_steady_HEADERS} DESTINATION ${ROM_INCDIR}/lspg_steady)
INSTALL(FILES ${rom_lspg_steady_generator_HEADERS} DESTINATION ${ROM_INCDIR}/lspg_steady/generator)


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
