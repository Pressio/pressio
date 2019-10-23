
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "model_for_static_checks_for_implicit_stepper.hpp"

TEST(ode_implicit, checkStepSetter){
  using namespace pressio;

  using sc_t = double;

  const auto dt_setter = [](const int & step, const sc_t & time, sc_t & dt) -> void{
			   dt = 1.1234;
			 };
  using setter_t = decltype(dt_setter);
  static_assert(ode::meta::is_legitimate_time_step_setter<setter_t, int, sc_t>::value, "");
}
