
execute_process(
  COMMAND ${CMAKE_COMMAND} --build . --target ${targetName}
  --config $<CONFIGURATION> ERROR_FILE ${targetLogFile})

#execute_process(COMMAND bash -c "")

execute_process(COMMAND bash -c "grep -r '${stringToGrepOne}' ${targetLogFile}" RESULT_VARIABLE Rone)
execute_process(COMMAND bash -c "grep -r '${stringToGrepTwo}' ${targetLogFile}" RESULT_VARIABLE Rtwo)
message("\n")
if(Rone)
  message(FATAL_ERROR "Expected string: ${stringToGrepOne}, not found")
elseif(Rtwo)
  message(FATAL_ERROR "Expected string: ${stringToGrepTwo}, not found")
else()
  message("This test is testing that an error occurs")
  message("I found the expected string: ${stringToGrepOne} and ${stringToGrepTwo}" )
endif()
