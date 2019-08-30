
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "rom_config.h"

## Set the base directory in order to find pressio_version.h
#ASSERT_DEFINED(PROJECT_SOURCE_DIR)
#INCLUDE_DIRECTORIES(${${PROJECT_NAME}_BINARY_DIR})

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

FILE(GLOB rom_gen_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

FILE(GLOB rom_experimental_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.hpp")
FILE(GLOB rom_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")

FILE(GLOB rom_decoder_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/decoder/*.hpp")
FILE(GLOB rom_decorators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/decorators/*.hpp")
FILE(GLOB rom_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/policies/*.hpp")

FILE(GLOB rom_galerkin_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/galerkin/*.hpp")
FILE(GLOB rom_galerkin_generator_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/galerkin/generator/*.hpp")

FILE(GLOB rom_lspg_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/*.hpp")
FILE(GLOB rom_lspg_unsteady_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/unsteady/*.hpp")
FILE(GLOB rom_lspg_unsteady_generator_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/unsteady/generator/*.hpp")
FILE(GLOB rom_lspg_steady_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/steady/*.hpp")
FILE(GLOB rom_lspg_steady_generator_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/lspg/steady/generator/*.hpp")

# sources
FILE(GLOB rom_gen_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
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
INSTALL(FILES ${rom_meta_HEADERS} DESTINATION ${ROM_INCDIR}/meta)

INSTALL(FILES ${rom_experimental_HEADERS} DESTINATION ${ROM_INCDIR}/experimental)

INSTALL(FILES ${rom_decoder_HEADERS} DESTINATION ${ROM_INCDIR}/decoder)
INSTALL(FILES ${rom_decorators_HEADERS} DESTINATION ${ROM_INCDIR}/decorators)
INSTALL(FILES ${rom_policies_HEADERS} DESTINATION ${ROM_INCDIR}/policies)

INSTALL(FILES ${rom_galerkin_HEADERS} DESTINATION ${ROM_INCDIR}/galerkin)
INSTALL(FILES ${rom_galerkin_generator_HEADERS} DESTINATION ${ROM_INCDIR}/galerkin/generator)

INSTALL(FILES ${rom_lspg_HEADERS} DESTINATION ${ROM_INCDIR}/lspg)
INSTALL(FILES ${rom_lspg_unsteady_HEADERS} DESTINATION ${ROM_INCDIR}/lspg/unsteady)
INSTALL(FILES ${rom_lspg_unsteady_generator_HEADERS} DESTINATION ${ROM_INCDIR}/lspg/unsteady/generator)
INSTALL(FILES ${rom_lspg_steady_HEADERS} DESTINATION ${ROM_INCDIR}/lspg/steady)
INSTALL(FILES ${rom_lspg_steady_generator_HEADERS} DESTINATION ${ROM_INCDIR}/lspg/steady/generator)


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
