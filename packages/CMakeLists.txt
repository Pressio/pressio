
# the directory where headers will be installed
set(INC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)

##################################################
# first, make list of all files in this directory
# including all the global includes files
##################################################

file(GLOB all_files_here "*")

foreach(f ${all_files_here})
  if(NOT f MATCHES "\\.txt"
  	 AND NOT f MATCHES ".DS_Store"
     AND NOT f MATCHES ".h.in"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/mpl"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/utils"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/containers"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/apps"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/qr"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/svd"
  	 AND NOT f MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/optimization"
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
install(DIRECTORY apps 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY qr 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY svd 		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY optimization	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY solvers 	DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ode		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY rom		DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
