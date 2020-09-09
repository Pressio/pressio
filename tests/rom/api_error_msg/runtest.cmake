
execute_process(
  COMMAND ${CMAKE_COMMAND} --build . --target ${targetName}
  --config $<CONFIGURATION> ERROR_FILE ${targetLogFile})

#execute_process(COMMAND bash -c "")

execute_process(COMMAND bash -c "cat ${targetLogFile}")
execute_process(COMMAND bash -c "grep -r '${stringToGrep}' ${targetLogFile}" RESULT_VARIABLE CMD_RESULT)
message("\n")
if(CMD_RESULT)
  message(FATAL_ERROR "Expected error message: ${stringToGrep}, not found")
else()
  message("This test is testing that an error occurs")
  message("I found the expected error message: ${stringToGrep}")
endif()
