
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"

namespace
{
struct FakeMassMatrixTypeForTesting{};
struct FakeStateTypeForTesting{};
struct FakeJacTypeForTesting{};

//
// velo only
//
struct ValidSystemWithVelocity1{
  using time_type = float;
  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;

  state_type createState() const{ return state_type(); }
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity1{
  using time_type = float;
  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;

  // comment out createVelocity so that it is invalid
  // velocity_type createVelocity() const{ return velocity_type(); }

  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity2{
  // comment out time_type so that it is valid
  // using time_type = float;

  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

//
// velo and mass matrix
//
struct ValidSystemWithVelocityAndMassMatrix1{
  using time_type = float;
  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}

  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }
  void massMatrix(const state_type &, double time, mass_matrix_type &) const{}
};

//
// velo and jacobian
//
struct ValidSystemWithVelocityAndJacobian1{
  using time_type = float;
  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;
  using jacobian_type = FakeJacTypeForTesting;

  state_type createState() const{ return state_type(); }
  velocity_type createVelocity() const{ return velocity_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
};

//
// velo and jacobian and mass matrix
//
struct ValidSystemWithVelocityAndJacobianAndMassMatrix1{
  using time_type = float;
  using state_type = FakeStateTypeForTesting;
  using velocity_type = state_type;
  using jacobian_type = FakeJacTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  velocity_type createVelocity() const{ return velocity_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }
  void massMatrix(const state_type &, double time, mass_matrix_type &) const{}
};

}//end anonym namespace

TEST(ode, concepts_continuous_time_system)
{
  using namespace pressio::ode;

  static_assert
    (continuous_time_system_with_at_least_velocity<
     ValidSystemWithVelocity1>::value, "");
  static_assert
    (continuous_time_system_with_mass_matrix_and_at_least_velocity<
     ValidSystemWithVelocityAndMassMatrix1>::value, "");
  static_assert
    (continuous_time_system_with_user_provided_jacobian<
     ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert
    (continuous_time_system_with_mass_matrix_and_user_provided_jacobian<
     ValidSystemWithVelocityAndJacobianAndMassMatrix1>::value, "");


  static_assert
    (continuous_time_system_with_at_least_velocity<
     ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert
    (continuous_time_system_with_at_least_velocity<
     ValidSystemWithVelocityAndMassMatrix1>::value, "");

  static_assert
    (!continuous_time_system_with_mass_matrix_and_user_provided_jacobian<
     ValidSystemWithVelocity1>::value, "");

  static_assert
    (!continuous_time_system_with_mass_matrix_and_at_least_velocity<
     ValidSystemWithVelocity1>::value, "");
  static_assert
    (!continuous_time_system_with_mass_matrix_and_at_least_velocity<
     ValidSystemWithVelocityAndJacobian1>::value, "");

  static_assert
    (!continuous_time_system_with_at_least_velocity<
     InvalidSystemWithVelocity1>::value, "");
  static_assert
    (!continuous_time_system_with_at_least_velocity<
     InvalidSystemWithVelocity2>::value, "");

  static_assert
    (!continuous_time_system_with_user_provided_jacobian<
     ValidSystemWithVelocity1>::value, "");
}
