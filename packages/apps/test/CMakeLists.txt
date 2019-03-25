
#------------------------------------------------------------------------
# to include the lumped headers to have access to package classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

ADD_SUBDIRECTORIES(burgers1d)
ADD_SUBDIRECTORIES(steady_adv_diff_1d)
