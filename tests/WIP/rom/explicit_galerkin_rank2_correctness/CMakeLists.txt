
set(name forward_euler_rank2_galerkin_correctness)
add_serial_exe_and_test(${name} rom main.cc "PASSED")
target_compile_definitions(rom_${name} PRIVATE EULER)

set(name rk4_rank2_galerkin_correctness)
add_serial_exe_and_test(${name} rom main.cc "PASSED")
target_compile_definitions(rom_${name} PRIVATE RK4)
