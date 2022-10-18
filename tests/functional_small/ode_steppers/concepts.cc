
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

struct FakeIndVarTypeForTesting{
  operator double(){ return double{}; }
};

//
// rhs only
//
struct System1{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void operator()(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  right_hand_side_type &    /*unused*/) const{}
};

struct System2{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;

  // omit createRightHandSide so that it is invalid

  void operator()(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  right_hand_side_type &    /*unused*/) const{}
};

struct System3{
  // omit independent_variable_type so that it is valid
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }

  void operator()(const state_type &       /*unused*/,
		  FakeIndVarTypeForTesting /*unused*/,
		  right_hand_side_type &   /*unused*/) const{}
};

//
// rhs and mass matrix
//
struct System4{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;
  using mass_matrix_type = Eigen::SparseMatrix<double>;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void operator()(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  right_hand_side_type &    /*unused*/,
		  mass_matrix_type &        /*unused*/) const{}
};

//
// rhs and jacobian
//
struct System6{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }

  void operator()(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  right_hand_side_type &    /*unused*/,
		  jacobian_type &           /*unused*/,
		  bool computeJac) const{}
};

//
// rhs and jacobian and mass matrix
//
struct System7{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using right_hand_side_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
  using mass_matrix_type = Eigen::SparseMatrix<double>;

  state_type createState() const{ return state_type(); }
  right_hand_side_type createRightHandSide() const{ return right_hand_side_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void operator()(const state_type &        /*unused*/,
		  independent_variable_type /*unused*/,
		  right_hand_side_type &    /*unused*/,
		  mass_matrix_type &        /*unused*/,
		  jacobian_type &           /*unused*/,
		  bool computeJac) const{}
};

//
// fully discrete
//
struct System9{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using discrete_residual_type = Eigen::VectorXd;
  using discrete_jacobian_type = Eigen::MatrixXd;

  state_type createState() const{ return state_type(); }

  discrete_residual_type createDiscreteResidual() const{
    return discrete_residual_type(); }
  discrete_jacobian_type createDiscreteJacobian() const{
    return discrete_jacobian_type(); }

  template<class StepCountType>
  void discreteResidualAndJacobian(StepCountType,
				   independent_variable_type /*unused*/,
				   independent_variable_type /*unused*/,
				   discrete_residual_type &  /*unused*/,
				   discrete_jacobian_type &  /*unused*/,
				   bool computeJacobian,
				   const state_type &        /*unused*/) const{}

  template<class StepCountType>
  void discreteResidualAndJacobian(StepCountType,
				   independent_variable_type /*unused*/,
				   independent_variable_type /*unused*/,
				   discrete_residual_type &  /*unused*/,
				   discrete_jacobian_type &  /*unused*/,
				   bool computeJacobian,
				   const state_type &        /*unused*/,
				   const state_type &        /*unused*/) const{}
};

// template<typename StateType, typename ResidualType>
// class ResidualPolicy1
// {
// public:
//   using independent_variable_type = double;
//   using state_type = StateType;
//   using residual_type = ResidualType;

//   state_type createState() const{ return state_type(); }
//   residual_type createResidual() const{ return residual_type(); }

//   template <typename prev_states_type, class rhs_container>
//   void operator()(pressio::ode::StepScheme          /*unused*/,
// 		  const StateType &		    /*unused*/,
// 		  const prev_states_type &	    /*unused*/,
// 		  rhs_container &		    /*unused*/,
// 		  const ::pressio::ode::StepEndAt<double> & /*unused*/,
// 		  ::pressio::ode::StepCount /*unused*/,
// 		  const ::pressio::ode::StepSize<double> & /*unused*/,
// 		  residual_type &		    /*unused*/) const
//   {}
// };

// template<typename StateType, typename JacobianType>
// class JacobianPolicy1
// {
// public:
//   using independent_variable_type = double;
//   using state_type = StateType;
//   using jacobian_type = JacobianType;

//   state_type createState() const{ return state_type(); }
//   jacobian_type createJacobian()   const{ return jacobian_type(); }

//   template <typename prev_states_type>
//   void operator()(pressio::ode::StepScheme          /*unused*/,
// 		  const StateType &		    /*unused*/,
// 		  const prev_states_type &	    /*unused*/,
// 		  const ::pressio::ode::StepEndAt<double> & /*unused*/,
// 		  ::pressio::ode::StepCount /*unused*/,
// 		  const ::pressio::ode::StepSize<double> & /*unused*/,
// 		  jacobian_type &		    /*unused*/) const
//   {}
// };


TEST(ode, concepts)
{
  using namespace pressio::ode;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert( SystemWithRhs<System1>, "");
  static_assert(!SystemWithRhsAndJacobian<System1>, "");

  static_assert(!SystemWithRhs<System2>, "");
  static_assert(!SystemWithRhsAndJacobian<System2>, "");

  static_assert(!SystemWithRhs<System3>, "");
  static_assert(!SystemWithRhsAndJacobian<System3>, "");

  static_assert( SystemWithRhsAndMassMatrix<System4>, "");

  static_assert( SystemWithRhsAndJacobian<System6>, "");
  static_assert( SystemWithRhsJacobianMassMatrix<System7>, "");

  static_assert(FullyDiscreteSystemWithJacobian<System9,  1>, "");
  static_assert(FullyDiscreteSystemWithJacobian<System9,  2>, "");
  static_assert(!FullyDiscreteSystemWithJacobian<System9, 3>, "");
  static_assert(!SystemWithRhs<System9>, "");
  static_assert(!SystemWithRhsAndJacobian<System9>, "");

#else

  static_assert(SystemWithRhs<System1>::value, "");
  static_assert(!SystemWithRhsAndJacobian<System1>::value, "");
  static_assert(!SystemWithRhs<System2>::value, "");
  static_assert(!SystemWithRhsAndJacobian<System2>::value, "");
  static_assert(!SystemWithRhs<System3>::value, "");
  static_assert(!SystemWithRhsAndJacobian<System3>::value, "");
  static_assert(SystemWithRhsAndMassMatrix<System4>::value, "");
  static_assert(SystemWithRhsAndJacobian<System6>::value, "");
  static_assert(SystemWithRhsJacobianMassMatrix<System7>::value, "");

  static_assert(FullyDiscreteSystemWithJacobian<System9,  1>::value, "");
  static_assert(FullyDiscreteSystemWithJacobian<System9,  2>::value, "");
  static_assert(!FullyDiscreteSystemWithJacobian<System9, 3>::value, "");
  static_assert(!SystemWithRhs<System9>::value, "");
  static_assert(!SystemWithRhsAndJacobian<System9>::value, "");
#endif

}
