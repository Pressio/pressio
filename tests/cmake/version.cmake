
# versioning
#=====================================================================
# adapted from Eigen
file(READ "${CMAKE_SOURCE_DIR}/include/pressio/pressio_macros.hpp" _pressio_macros)

string(REGEX MATCH "define[ \t]+PRESSIO_MAJOR_VERSION[ \t]+([0-9]+)" _pressio_major_version_match "${_pressio_macros}")
set(PRESSIO_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+PRESSIO_MINOR_VERSION[ \t]+([0-9]+)" _pressio_minor_version_match "${_pressio_macros}")
set(PRESSIO_MINOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+PRESSIO_PATCH_VERSION[ \t]+([0-9]+)" _pressio_patch_version_match "${_pressio_macros}")
set(PRESSIO_PATCH_VERSION "${CMAKE_MATCH_1}")
set(PRESSIO_VERSION_NUMBER ${PRESSIO_MAJOR_VERSION}.${PRESSIO_MINOR_VERSION}.${PRESSIO_PATCH_VERSION})
message("${Magenta}>> PRESSIO: version = ${PRESSIO_VERSION_NUMBER} ${ColourReset}")
