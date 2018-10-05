
FILE(GLOB tinympl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

FILE(GLOB variadic_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/variadic/*.hpp")

install (FILES ${tinympl_HEADERS} DESTINATION include/tinympl)
install (FILES ${variadic_HEADERS} DESTINATION include/tinympl/variadic)
