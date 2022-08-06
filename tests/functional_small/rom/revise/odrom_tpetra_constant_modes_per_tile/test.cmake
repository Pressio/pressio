include(FindUnixCommands)

set(CMD "mpirun --oversubscribe -np 5 ${EXENAME}")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "executable failed")
else()
  message("exe succeeded!")
endif()

set(CMD "python3 check_rec.py")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "reconstruction failed")
else()
  message("reconstruction succeeded!")
endif()

set(CMD "python3 check_proj.py")
execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "projection failed")
else()
  message("projection succeeded!")
endif()