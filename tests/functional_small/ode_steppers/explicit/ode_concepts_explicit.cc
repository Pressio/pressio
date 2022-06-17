
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
  using independent_variable_type = float;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}
};

struct InvalidSystemWithVelocity1{
  using independent_variable_type = float;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;

  // comment out createRightHandSide so that it is invalid
  // right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}
};

struct InvalidSystemWithVelocity2{
  // comment out independent_variable_type so that it is valid
  // using independent_variable_type = float;

  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}
};

//
// velo and mass matrix
//
struct ValidSystemWithVelocityAndMassMatrix1{
  using independent_variable_type = float;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}

  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }
  void massMatrix(const state_type &, double time, mass_matrix_type &) const{}
};

//
// velo and jacobian
//
struct ValidSystemWithVelocityAndJacobian1{
  using independent_variable_type = float;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;
  using jacobian_type = FakeJacTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
};

//
// velo and jacobian and mass matrix
//
struct ValidSystemWithVelocityAndJacobianAndMassMatrix1{
  using independent_variable_type = float;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = state_type;
  using jacobian_type = FakeJacTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void rightHandSide(const state_type &, double time, right_hand_side_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }
  void massMatrix(const state_type &, double time, mass_matrix_type &) const{}
};

}//end anonym namespace

TEST(ode, concepts_continuous_time_system)
{
  using namespace pressio::ode;

  static_assert
    (SemiDiscreteSystem<
     ValidSystemWithVelocity1>::value, "");
  static_assert
    (SemiDiscreteSystemWithMassMatrix<
     ValidSystemWithVelocityAndMassMatrix1>::value, "");
  static_assert
    (SemiDiscreteSystemWithJacobian<
     ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert
    (SemiDiscreteSystemComplete<
     ValidSystemWithVelocityAndJacobianAndMassMatrix1>::value, "");


  static_assert
    (SemiDiscreteSystem<
     ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert
    (SemiDiscreteSystem<
     ValidSystemWithVelocityAndMassMatrix1>::value, "");

  static_assert
    (!SemiDiscreteSystemComplete<
     ValidSystemWithVelocity1>::value, "");

  static_assert
    (!SemiDiscreteSystemWithMassMatrix<
     ValidSystemWithVelocity1>::value, "");
  static_assert
    (!SemiDiscreteSystemWithMassMatrix<
     ValidSystemWithVelocityAndJacobian1>::value, "");

  static_assert
    (!SemiDiscreteSystem<
     InvalidSystemWithVelocity1>::value, "");
  static_assert
    (!SemiDiscreteSystem<
     InvalidSystemWithVelocity2>::value, "");

  static_assert
    (!SemiDiscreteSystemWithJacobian<
     ValidSystemWithVelocity1>::value, "");

}
