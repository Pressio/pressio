
#include <gtest/gtest.h>

#include "pressio/type_traits.hpp"

struct FakeStateTypeForTesting{};
struct FakeRhsTypeForTesting{};
struct FakeMassMatrixTypeForTesting{};
struct FakeJacTypeForTesting{};
struct FakeIndVarTypeForTesting{
  operator double(){ return double{}; }
};

namespace pressio{
template<> struct Traits<FakeStateTypeForTesting>{ using scalar_type = double; };
template<> struct Traits<FakeRhsTypeForTesting>{ using scalar_type = double; };
template<> struct Traits<FakeMassMatrixTypeForTesting>{ using scalar_type = double; };
template<> struct Traits<FakeJacTypeForTesting>{ using scalar_type = double; };
}

#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_steppers_implicit.hpp"

//
// rhs only
//
struct System1{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}
};

struct System2{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;

  // omit createRightHandSide so that it is invalid

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}
};

struct System3{
  // omit independent_variable_type so that it is valid
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void rightHandSide(const state_type &       /*unused*/,
		     FakeIndVarTypeForTesting /*unused*/,
		     right_hand_side_type &   /*unused*/) const{}
};

//
// rhs and mass matrix
//
struct System4{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}

  void massMatrix(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  mass_matrix_type &        /*unused*/) const{}
};

//
// rhs and constant mass matrix
//
struct System5{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}

  void massMatrix(mass_matrix_type & /*unused*/) const{}
};

//
// rhs and jacobian
//
struct System6{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  using jacobian_type = FakeJacTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}

  void jacobian(const state_type &        /*unused*/,
		independent_variable_type /*unused*/,
		jacobian_type &           /*unused*/) const{}
};

//
// rhs and jacobian and mass matrix
//
struct System7{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  using jacobian_type = FakeJacTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}

  void massMatrix(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  mass_matrix_type &        /*unused*/) const{}

  void jacobian(const state_type &        /*unused*/,
		independent_variable_type /*unused*/,
		jacobian_type &           /*unused*/) const{}
};

//
// rhs and jacobian and constant mass matrix
//
struct System8{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using right_hand_side_type = FakeRhsTypeForTesting;
  using jacobian_type = FakeJacTypeForTesting;
  using mass_matrix_type = FakeMassMatrixTypeForTesting;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void rightHandSide(const state_type &        /*unused*/,
		     independent_variable_type /*unused*/,
		     right_hand_side_type &    /*unused*/) const{}

  void massMatrix(mass_matrix_type &        /*unused*/) const{}

  void jacobian(const state_type &        /*unused*/,
		independent_variable_type /*unused*/,
		jacobian_type &           /*unused*/) const{}
};


struct System9{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = FakeStateTypeForTesting;
  using discrete_residual_type = FakeRhsTypeForTesting;
  using discrete_jacobian_type = FakeJacTypeForTesting;

  state_type createState() const{ return state_type(); }

  discrete_residual_type createDiscreteResidual() const{
    return discrete_residual_type(); }
  discrete_jacobian_type createDiscreteJacobian() const{
    return discrete_jacobian_type(); }

  template<class StepCountType>
  void discreteResidual(StepCountType,
			independent_variable_type /*unused*/,
			independent_variable_type /*unused*/,
			discrete_residual_type &  /*unused*/,
			const state_type &        /*unused*/) const{}

  template<class StepCountType>
  void discreteResidual(StepCountType,
			independent_variable_type /*unused*/,
			independent_variable_type /*unused*/,
			discrete_residual_type &  /*unused*/,
			const state_type &        /*unused*/,
			const state_type &        /*unused*/) const{}

  template<class StepCountType>
  void discreteJacobian(StepCountType,
			independent_variable_type /*unused*/,
			independent_variable_type /*unused*/,
			discrete_jacobian_type &  /*unused*/,
			const state_type &        /*unused*/) const{}

  template<class StepCountType>
  void discreteJacobian(StepCountType,
			independent_variable_type /*unused*/,
			independent_variable_type /*unused*/,
			discrete_jacobian_type &  /*unused*/,
			const state_type &        /*unused*/,
			const state_type &        /*unused*/) const{}
};

template<typename StateType, typename ResidualType>
class ResidualPolicy1
{
public:
  using independent_variable_type = double;
  using state_type = StateType;
  using residual_type = ResidualType;

  state_type createState() const{ return state_type(); }
  residual_type createResidual() const{ return residual_type(); }

  template <typename prev_states_type, class rhs_container>
  void operator()(pressio::ode::StepScheme          /*unused*/,
		  const StateType &		    /*unused*/,
		  const prev_states_type &	    /*unused*/,
		  rhs_container &		    /*unused*/,
		  const ::pressio::ode::StepEndAt<double> & /*unused*/,
		  ::pressio::ode::StepCount /*unused*/,
		  const ::pressio::ode::StepSize<double> & /*unused*/,
		  residual_type &		    /*unused*/) const
  {}
};

template<typename StateType, typename JacobianType>
class JacobianPolicy1
{
public:
  using independent_variable_type = double;
  using state_type = StateType;
  using jacobian_type = JacobianType;

  state_type createState() const{ return state_type(); }
  jacobian_type createJacobian()   const{ return jacobian_type(); }

  template <typename prev_states_type>
  void operator()(pressio::ode::StepScheme          /*unused*/,
		  const StateType &		    /*unused*/,
		  const prev_states_type &	    /*unused*/,
		  const ::pressio::ode::StepEndAt<double> & /*unused*/,
		  ::pressio::ode::StepCount /*unused*/,
		  const ::pressio::ode::StepSize<double> & /*unused*/,
		  jacobian_type &		    /*unused*/) const
  {}
};


TEST(ode, concepts)
{
  using namespace pressio::ode;

  static_assert( OdeSystem<System1>::value, "");
  static_assert(!OdeSystemWithJacobian<System1>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System1>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System1>::value, "");
  static_assert(!OdeSystemComplete<System1>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System1>::value, "");

  static_assert(!OdeSystem<System2>::value, "");
  static_assert(!OdeSystemWithJacobian<System2>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System2>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System2>::value, "");
  static_assert(!OdeSystemComplete<System2>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System2>::value, "");

  static_assert(!OdeSystem<System3>::value, "");
  static_assert(!OdeSystemWithJacobian<System3>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System3>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System3>::value, "");
  static_assert(!OdeSystemComplete<System3>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System3>::value, "");

  static_assert( OdeSystem<System4>::value, "");
  static_assert(!OdeSystemWithJacobian<System4>::value, "");
  static_assert( OdeSystemWithMassMatrix<System4>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System4>::value, "");
  static_assert(!OdeSystemComplete<System4>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System4>::value, "");

  static_assert( OdeSystem<System5>::value, "");
  static_assert(!OdeSystemWithJacobian<System5>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System5>::value, "");
  static_assert( OdeSystemWithConstantMassMatrix<System5>::value, "");
  static_assert(!OdeSystemComplete<System5>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System5>::value, "");

  static_assert( OdeSystem<System6>::value, "");
  static_assert( OdeSystemWithJacobian<System6>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System6>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System6>::value, "");
  static_assert(!OdeSystemComplete<System6>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System6>::value, "");

  static_assert( OdeSystem<System7>::value, "");
  static_assert( OdeSystemWithJacobian<System7>::value, "");
  static_assert( OdeSystemWithMassMatrix<System7>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System7>::value, "");
  static_assert( OdeSystemComplete<System7>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System7>::value, "");

  static_assert( OdeSystem<System8>::value, "");
  static_assert( OdeSystemWithJacobian<System8>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System8>::value, "");
  static_assert( OdeSystemWithConstantMassMatrix<System8>::value, "");
  static_assert(!OdeSystemComplete<System8>::value, "");
  static_assert( OdeSystemCompleteWithConstantMassMatrix<System8>::value, "");

  static_assert(FullyDiscreteSystemWithJacobian<System9,  1>::value, "");
  static_assert(FullyDiscreteSystemWithJacobian<System9,  2>::value, "");
  static_assert(!FullyDiscreteSystemWithJacobian<System9, 3>::value, "");
  static_assert(!OdeSystem<System9>::value, "");
  static_assert(!OdeSystemWithJacobian<System9>::value, "");
  static_assert(!OdeSystemWithMassMatrix<System9>::value, "");
  static_assert(!OdeSystemWithConstantMassMatrix<System9>::value, "");
  static_assert(!OdeSystemComplete<System9>::value, "");
  static_assert(!OdeSystemCompleteWithConstantMassMatrix<System9>::value, "");

  {
    using state_t = FakeStateTypeForTesting;
    using res_t   = FakeStateTypeForTesting;
    using jac_t   = FakeJacTypeForTesting;
    static_assert(ImplicitEulerResidualPolicy<ResidualPolicy1<state_t, res_t>>::value, "");
    static_assert(ImplicitEulerJacobianPolicy<JacobianPolicy1<state_t, jac_t>>::value, "");
  }
}
