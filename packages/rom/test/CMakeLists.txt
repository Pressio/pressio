
#------------------------------------------------------------------------
# to include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

# include the utilities in utils (adding this here means
# all children cmakelists will inherit this include )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utils)

ADD_SUBDIRECTORIES(burgers1d steady_linear_adv_diff_2d)
