
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "model_for_static_checks_for_implicit_stepper.hpp"

TEST(ode_implicit, checkModelsScalar){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingScalarTypedef;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::containers::meta::has_scalar_typedef<app_t>::value, "");
}

TEST(ode_implicit, checkModelsStateMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingStateTypedef;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::has_state_typedef<app_t>::value, "");
}

TEST(ode_implicit, checkModelsVelocityMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingVelocityTypedef;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::has_velocity_typedef<app_t>::value, "");
}


TEST(ode_implicit, checkModelsJacobianMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingJacobianTypedef;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::has_jacobian_typedef<app_t>::value, "");
}

TEST(ode_implicit, checkModelsVelocityMethodsMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingVelocityMethods;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using nvel_t   = typename app_t::velocity_type;

  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_velocity_methods<
		app_t, nstate_t, nvel_t, sc_t>::value, "");
}

TEST(ode_implicit, checkModelsVelocityNonVoidMethodMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingVelocityNonVoidMethod;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using nvel_t   = typename app_t::velocity_type;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_velocity_methods<
		app_t, nstate_t, nvel_t, sc_t>::value, "");
}

TEST(ode_implicit, checkModelsVelocityVoidMethodMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingVelocityVoidMethod;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using nvel_t   = typename app_t::velocity_type;
  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_velocity_methods<
		app_t, nstate_t, nvel_t, sc_t>::value, "");
}


TEST(ode_implicit, checkModelsJacobianMethodsMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingJacobianMethods;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using njac_t   = typename app_t::jacobian_type;

  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_jacobian_methods<
		app_t, nstate_t, njac_t, sc_t>::value, "");
}

TEST(ode_implicit, checkModelsJacobianVoidMethodMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingJacobianVoidMethod;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using njac_t   = typename app_t::jacobian_type;

  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_jacobian_methods<
		app_t, nstate_t, njac_t, sc_t>::value, "");
}

TEST(ode_implicit, checkModelsJacobianNonVoidMethodMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingJacobianNonVoidMethod;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using njac_t   = typename app_t::jacobian_type;

  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
  static_assert(!::pressio::ode::meta::model_has_needed_jacobian_methods<
		app_t, nstate_t, njac_t, sc_t>::value, "");
}

TEST(ode_implicit, checkModelsConstMissing){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitMissingConst;
  using sc_t	 = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  using njac_t   = typename app_t::jacobian_type;

  static_assert(!ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
}


TEST(ode_implicit, checkModelsValid){
  using namespace pressio;
  using app_t	 = ode::testing::ModelForImplicitValid;
  static_assert(ode::meta::is_legitimate_model_for_implicit_ode<app_t>::value, "");
}
