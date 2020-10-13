
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"

TEST(user_defined_model, admissibleExplicitOde){
  using namespace pressio;

  // struct fakeapp{
  //   using scalar_type = double;
  //   using state_type = std::vector<double>;
  //   using residual_type = std::vector<double>;

  //   void residual(const state_type & y,
  //   		  space_residual_type & R,
  //   		  double t){
  //   };

  //   state_type residual(const state_type & y,
  // 			scalar_type y2){
  //     state_type a;
  //     return a;
  //   };
  // };

  using app_t = ::pressio::ode::testing::fakeAppForTraitsForExp;
  static_assert(
    ode::concepts::continuous_time_system_without_jacobian<app_t>::value, "");
}
