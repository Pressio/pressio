
#----------------------------------------------------------------
# A) Define the package
#----------------------------------------------------------------
TRIBITS_PACKAGE(rom ENABLE_SHADOWING_WARNINGS CLEANED )

#----------------------------------------------------------------
# B) Platform-specific checks
#----------------------------------------------------------------
INCLUDE(CheckFor__int64)
#CHECK_FOR___INT64(HAVE_SIMPLECXX___INT64)

#----------------------------------------------------------------
# C) Set up package-specific options
#----------------------------------------------------------------
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#----------------------------------------------------------------
# D) Add the libraries, tests, and examples
#----------------------------------------------------------------
# needs headers from core, solvers, ode, apps, optimization 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../core/lumped_include_headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../solvers/lumped_include_headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../svd/lumped_include_headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../optimization/lumped_include_headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../ode/lumped_include_headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../apps/lumped_include_headers)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unit_tests)
TRIBITS_ADD_TEST_DIRECTORIES(test)
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#----------------------------------------------------------------
# E) Do standard post processing
#----------------------------------------------------------------
TRIBITS_PACKAGE_POSTPROCESS()
