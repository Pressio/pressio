
#------------------------------------------------------------------------
# to include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

add_serial_unit_test(meta_basic)
add_serial_unit_test(meta_detect_operators)
