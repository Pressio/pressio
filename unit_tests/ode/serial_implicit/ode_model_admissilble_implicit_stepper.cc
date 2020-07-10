
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "model_for_static_checks_for_implicit_stepper.hpp"

TEST(ode_implicit, checkModelsScalar){
  using namespace pressio;
  using system_t	 = ode::testing::ModelForImplicitMissingScalarTypedef;
  static_assert(!ode::concepts::continuous_time_implicit_system<system_t>::value, "");
  static_assert(!ode::concepts::discrete_time_system_implicit_stepping<system_t>::value, "");
  static_assert(!::pressio::containers::predicates::has_scalar_typedef<system_t>::value, "");
}

TEST(ode_implicit, checkModelsStateMissing){
  using namespace pressio;
  using system_t	 = ode::testing::ModelForImplicitMissingStateTypedef;
  static_assert(!ode::concepts::continuous_time_implicit_system<system_t>::value, "");
  static_assert(!ode::concepts::discrete_time_system_implicit_stepping<system_t>::value, "");
  static_assert(!::pressio::ode::predicates::has_state_typedef<system_t>::value, "");
}

TEST(ode_implicit, checkModelsVelocityMissing){
  using namespace pressio;
  using system_t	 = ode::testing::ModelForImplicitMissingVelocityTypedef;
  static_assert(!ode::concepts::continuous_time_implicit_system<system_t>::value, "");
  static_assert(!ode::concepts::discrete_time_system_implicit_stepping<system_t>::value, "");
  static_assert(!::pressio::ode::predicates::has_velocity_typedef<system_t>::value, "");
}

TEST(ode_implicit, checkModelsJacobianMissing){
  using namespace pressio;
  using system_t	 = ode::testing::ModelForImplicitMissingJacobianTypedef;
  static_assert(!ode::concepts::continuous_time_implicit_system<system_t>::value, "");
  static_assert(!ode::concepts::discrete_time_system_implicit_stepping<system_t>::value, "");
  static_assert(!::pressio::ode::predicates::has_jacobian_typedef<system_t>::value, "");
}
