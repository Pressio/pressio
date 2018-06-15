
# Set the package name
TRIBITS_PACKAGE(PyROM DISABLE_STRONG_WARNINGS)

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "PyROM can only be built with shared libraries. Building of shared libraries is currently set to OFF. To enable shared libraries please set the cache variable \"BUILD_SHARED_LIBS\" to ON")
ENDIF()

# Set the package version number
# SET(PyROM_VERSION 12.13)

# Define the list of PyROM packages and modules.  Packages are
# the top-level names only, such as Epetra, AztecOO and NOX.  The list
# of modules is similar to the packages, except it contains the names
# of python modules that correspond to nested namespaces, such as
# NOX.__init__ or NOX.Epetra.Interface.
SET(PyROM_PACKAGES "")
SET(PyROM_MODULES  "")

APPEND_SET(PyROM_PACKAGES Teuchos    )
APPEND_SET(PyROM_MODULES  Teuchos    )

# IF(PyROM_ENABLE_Epetra)
#   APPEND_SET(PyROM_PACKAGES Epetra)
#   APPEND_SET(PyROM_MODULES  Epetra)
# ENDIF(PyROM_ENABLE_Epetra)

# IF(PyROM_ENABLE_Isorropia)
#   APPEND_SET(PyROM_PACKAGES Isorropia         )
#   APPEND_SET(PyROM_MODULES  Isorropia.__init__)
#   IF(PyROM_ENABLE_Epetra)
#     Append_SET(PyROM_MODULES Isorropia.Epetra)
#   ENDIF(PyROM_ENABLE_Epetra)
# ENDIF(PyROM_ENABLE_Isorropia)

# Python interpreter and header files are required
FIND_PACKAGE(PythonInterp 2.7 REQUIRED)
FIND_PACKAGE(PythonInclude REQUIRED)

# Python module numpy is required
FIND_PACKAGE(NumPy 1.0 REQUIRED)

# Simple wrapper interface generator (SWIG) is required. SWIG 3.0.11
# is required to support the %module(moduleimport=...) option
SET(SWIG_VERSION_FLOOR   3.0.11)
FIND_PACKAGE(SWIG REQUIRED)
IF(${SWIG_VERSION} VERSION_LESS ${SWIG_VERSION_FLOOR})
  MESSAGE(FATAL_ERROR
    "SWIG version " ${SWIG_VERSION}
    " is less than minimum version " ${SWIG_VERSION_FLOOR}
    )
ENDIF()
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)

# Python module mpi4py is optional, and only of interest if MPI is
# enabled
IF(TPL_ENABLE_MPI)
  SET(TPL_ENABLE_MPI4PY ON CACHE BOOL "Use Mpi4Py")
  IF(TPL_ENABLE_MPI4PY)
    FIND_PACKAGE(Mpi4Py)
  ELSE(TPL_ENABLE_MPI4PY)
    SET(Mpi4Py_FOUND OFF)
  ENDIF(TPL_ENABLE_MPI4PY)
ENDIF(TPL_ENABLE_MPI)

# # Python module petsc4py is required if NOX_ENABLE_PETSC is true
# IF(NOX_ENABLE_PETSC)
#   FIND_PACKAGE(PETSc4Py REQUIRED)
# ENDIF(NOX_ENABLE_PETSC)

# Print the enabled PyROM packages
MESSAGE(STATUS "Enabled PyROM packages:")
MESSAGE(STATUS "    ${PyROM_PACKAGES}")

# Determine the PyROM default installation directory prefix.  If
# CMAKE_INSTALL_PREFIX is set by the user, then use it.  Else the
# prefix is set to the value of python's sys.prefix.
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.prefix)"
  OUTPUT_VARIABLE PYTHON_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyROM_DEFAULT_INSTALL_PREFIX ${PYTHON_PREFIX})
ELSE(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyROM_DEFAULT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Set the PyROM install prefix
SET(PyROM_INSTALL_PREFIX ${PyROM_DEFAULT_INSTALL_PREFIX}
  CACHE PATH "The path prefix for where PyROM will be installed, e.g. /usr/local")

# # The default status for whether we build docstrings for PyROM is
# # based on whether the user has doxygen installed.
# FIND_PACKAGE(Doxygen)
# IF(DOXYGEN_FOUND)
#   IF(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND PyROM_DOCSTRINGS)
#     MESSAGE(FATAL_ERROR 
#       "Doxygen version ${DOXYGEN_VERSION} has a bug that makes it incompatible with "
#       "PyROM. To work around this, you may disable PyROM docstrings by setting "
#       "PyROM_DOCSTRINGS=OFF. "
#       "Alternatively, you may choose to use a different version of doxygen. If you "
#       "choose this route, be sure that the new version of doxygen is found first in "
#       "your environment path.")
#   ELSE(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND PyROM_DOCSTRINGS)
#     SET(PyROM_DOCSTRINGS_DEFAULT ON)
#   ENDIF(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND PyROM_DOCSTRINGS)
# ELSE(DOXYGEN_FOUND)
#   SET(PyROM_DOCSTRINGS_DEFAULT OFF)
#   MESSAGE(STATUS "doxygen not found!")
# ENDIF(DOXYGEN_FOUND)

# # Set the PyROM docstrings flag
# SET(PyROM_DOCSTRINGS ${PyROM_DOCSTRINGS_DEFAULT} CACHE BOOL
#     "a flag that determines whether Python docstrings will be built using doxygen")
# IF(PyROM_DOCSTRINGS AND NOT DOXYGEN_FOUND)
#   MESSAGE(ERROR "PyROM docstrings enabled but doxygen not found")
# ENDIF(PyROM_DOCSTRINGS AND NOT DOXYGEN_FOUND)

# # Add the util directory
# ADD_SUBDIRECTORY(util)

# # Add the doc directory, which builds the python documentation strings
# ADD_SUBDIRECTORY(doc)

# Add the src directory, which builds the PyROM library,
# generates the wrapper code, builds the PyROM python package and
# compiles the PyROM extension modules
ADD_SUBDIRECTORY(src)

# Add the test scripts
TRIBITS_ADD_TEST_DIRECTORIES(test)

# Add the example scripts
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

# Execute the package postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
