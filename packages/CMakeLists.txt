
# the directory where headers will be installed
set(INC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)

##################################################
# first, make list of all files in this directory
# including all the global includes files
##################################################

# Dependency structure
# mpl/utils/containers/ops/qr/svd/optimizers/solvers/ode/rom
#
# the apps package is not a real dependency of rom if
# pressio needs to be used as a library because apps is
# only used for testing
#
# The dependency structure is reflected in the order below.
#

file(GLOB all_files_here "*")

foreach(f ${all_files_here})
  if(NOT f MATCHES "\\.txt"
      AND NOT f MATCHES ".DS_Store"
      AND NOT f MATCHES ".h.in"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/mpl"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/utils"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/containers"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/ops"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/apps"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/qr"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/svd"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/optimizers"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/solvers"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/ode"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/rom"
      AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/[.].+")
    list(APPEND files_to_install ${f})
  endif()
endforeach(f ${all_files_here})


##################################################
# install the packages headers
##################################################

install(FILES
  ${files_to_install}
  DESTINATION ${INC_INSTALL_DIR} #COMPONENT Devel
)

##################################################
# install the src code from individual packages
##################################################

install(DIRECTORY mpl		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY utils 	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY containers	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ops		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY apps 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY qr 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY svd 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY optimizers	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY solvers 	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ode		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY rom		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
