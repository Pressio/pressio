
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "ODE_ALL"
#include "reference_apps_for_testing.hpp"

TEST(user_defined_model, admissibleExplicitOde){
  using namespace rompp;

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

  using app_t = ::rompp::ode::testing::fakeAppForTraitsForExp;
  // static_assert(::rompp::mpl::is_detected<
		// ode::meta::has_scalar_typedef,
		// app_t>::value, " ");

  // static_assert(core::meta::is_detected<
  // 		ode::meta::has_residual_method_callable_with_two_args,
  // 		app_t, typename fakeapp::state_type,
  // 		typename fakeapp::state_type>::value, "");

  // static_assert(core::meta::is_detected<
  // 		ode::meta::has_residual_method_callable_with_three_args,
  // 		app_t,
  // 		typename fakeapp::state_type,
  // 		typename fakeapp::state_type,
  // 		double>::value, "");

  static_assert(
    ode::meta::is_legitimate_model_for_explicit_ode<app_t>::value, "");
  static_assert(
    !ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");


}
