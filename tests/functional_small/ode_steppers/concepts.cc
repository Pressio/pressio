
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

// //
// // rhs and mass matrix
// //
// struct System4{
//   using independent_variable_type = FakeIndVarTypeForTesting;
//   using state_type = FakeStateTypeForTesting;
//   using right_hand_side_type = FakeRhsTypeForTesting;
//   using mass_matrix_type = FakeMassMatrixTypeForTesting;

//   state_type createState() const{ return state_type(); }
//   right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
//   mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

//   void rightHandSide(const state_type &        /*unused*/,
// 		     independent_variable_type /*unused*/,
// 		     right_hand_side_type &    /*unused*/) const{}

//   void massMatrix(const state_type &        /*unused*/,
// 		  independent_variable_type /*unused*/,
// 		  mass_matrix_type &        /*unused*/) const{}
// };

// //
// // rhs and constant mass matrix
// //
// struct System5{
//   using independent_variable_type = FakeIndVarTypeForTesting;
//   using state_type = FakeStateTypeForTesting;
//   using right_hand_side_type = FakeRhsTypeForTesting;
//   using mass_matrix_type = FakeMassMatrixTypeForTesting;

//   state_type createState() const{ return state_type(); }
//   right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
//   mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

//   void rightHandSide(const state_type &        /*unused*/,
// 		     independent_variable_type /*unused*/,
// 		     right_hand_side_type &    /*unused*/) const{}

//   void massMatrix(mass_matrix_type & /*unused*/) const{}
// };

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

// //
// // rhs and jacobian and mass matrix
// //
// struct System7{
//   using independent_variable_type = FakeIndVarTypeForTesting;
//   using state_type = FakeStateTypeForTesting;
//   using right_hand_side_type = FakeRhsTypeForTesting;
//   using jacobian_type = FakeJacTypeForTesting;
//   using mass_matrix_type = FakeMassMatrixTypeForTesting;

//   state_type createState() const{ return state_type(); }
//   right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
//   jacobian_type createJacobian() const{ return jacobian_type(); }
//   mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

//   void rightHandSide(const state_type &        /*unused*/,
// 		     independent_variable_type /*unused*/,
// 		     right_hand_side_type &    /*unused*/) const{}

//   void massMatrix(const state_type &        /*unused*/,
// 		  independent_variable_type /*unused*/,
// 		  mass_matrix_type &        /*unused*/) const{}

//   void jacobian(const state_type &        /*unused*/,
// 		independent_variable_type /*unused*/,
// 		jacobian_type &           /*unused*/) const{}
// };

// //
// // rhs and jacobian and constant mass matrix
// //
// struct System8{
//   using independent_variable_type = FakeIndVarTypeForTesting;
//   using state_type = FakeStateTypeForTesting;
//   using right_hand_side_type = FakeRhsTypeForTesting;
//   using jacobian_type = FakeJacTypeForTesting;
//   using mass_matrix_type = FakeMassMatrixTypeForTesting;

//   state_type createState() const{ return state_type(); }
//   right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
//   jacobian_type createJacobian() const{ return jacobian_type(); }
//   mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

//   void rightHandSide(const state_type &        /*unused*/,
// 		     independent_variable_type /*unused*/,
// 		     right_hand_side_type &    /*unused*/) const{}

//   void massMatrix(mass_matrix_type &        /*unused*/) const{}

//   void jacobian(const state_type &        /*unused*/,
// 		independent_variable_type /*unused*/,
// 		jacobian_type &           /*unused*/) const{}
// };


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

  static_assert( OdeRhsEvaluator<System1>::value, "");
  static_assert(!OdeRhsAndJacobianEvaluator<System1>::value, "");

  static_assert(!OdeRhsEvaluator<System2>::value, "");
  static_assert(!OdeRhsAndJacobianEvaluator<System2>::value, "");

  static_assert(!OdeRhsEvaluator<System3>::value, "");
  static_assert(!OdeRhsAndJacobianEvaluator<System3>::value, "");

  static_assert( OdeRhsEvaluator<System6>::value, "");
  static_assert( OdeRhsAndJacobianEvaluator<System6>::value, "");

  static_assert(FullyDiscreteSystemWithJacobian<System9,  1>::value, "");
  static_assert(FullyDiscreteSystemWithJacobian<System9,  2>::value, "");
  static_assert(!FullyDiscreteSystemWithJacobian<System9, 3>::value, "");
  static_assert(!OdeRhsEvaluator<System9>::value, "");
  static_assert(!OdeRhsAndJacobianEvaluator<System9>::value, "");

  {
    using state_t = FakeStateTypeForTesting;
    using res_t   = FakeStateTypeForTesting;
    using jac_t   = FakeJacTypeForTesting;
    static_assert(ImplicitEulerResidualPolicy<ResidualPolicy1<state_t, res_t>>::value, "");
    static_assert(ImplicitEulerJacobianPolicy<JacobianPolicy1<state_t, jac_t>>::value, "");
  }
}
