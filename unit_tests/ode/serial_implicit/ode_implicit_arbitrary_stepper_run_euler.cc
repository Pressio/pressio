
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"

struct MyApp
{

  /*
    dy
    -- = -10*y
    dt

    y(0) = 1;
    y(1) = 2;
    y(2) = 3;
   */
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

private:
  template <typename step_t, typename state_type>
  void timeDiscreteResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    const auto f =  this->velocity(yn, time);
    R = yn - ynm1 - dt * f;
  }

  template <typename step_t, typename state_t>
  void timeDiscreteJacobianImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				jacobian_type & J,
				const state_t & yn,
				const state_t & ynm1) const
  {
    J =  this->jacobian(yn, time);
    constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    J.coeffs() *= -dt;
    J.coeffRef(0,0) += one;
    J.coeffRef(1,1) += one;
    J.coeffRef(2,2) += one;
  }

public:
  void velocity(const state_type & yIn,
		scalar_type t,
		velocity_type & f) const{
    f = -10. * yIn;
  };

  velocity_type velocity(const state_type & yIn,
			 scalar_type t) const{
    velocity_type f(3);
    this->velocity(yIn, t, f);
    return f;
  };

  void jacobian(const state_type & yIn,
  		scalar_type t,
		jacobian_type & JJ) const
  {
    typedef Eigen::Triplet<scalar_type> Tr;
    std::vector<Tr> tripletList;
    tripletList.push_back( Tr( 0, 0, -10.) );
    tripletList.push_back( Tr( 1, 1, -10.) );
    tripletList.push_back( Tr( 2, 2, -10.) );
    JJ.setFromTriplets(tripletList.begin(), tripletList.end());
  };

  jacobian_type jacobian(const state_type & yIn,
  			 scalar_type t) const{
    jacobian_type JJ(3,3);
    this->jacobian(yIn, t, JJ);
    return JJ;
  };
  //--------------------------------------------

public:
  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            residual_type & R,
                            Args && ... states) const
  {
    this->timeDiscreteResidualImpl( step, time, dt, R, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void timeDiscreteJacobian(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            jacobian_type & J,
                            Args && ... states) const
  {
    this->timeDiscreteJacobianImpl(step, time, dt, J, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  residual_type timeDiscreteResidual(const step_t & step,
                                     const scalar_type & time,
                                     const scalar_type & dt,
                                     Args && ... states) const
  {
    residual_type R(3);
    this->timeDiscreteResidual(step, time, dt, R, std::forward<Args>(states)...);
    return R;
  }

  template <typename step_t, typename ... Args>
  jacobian_type timeDiscreteJacobian(const step_t & step,
                                     const scalar_type & time,
                                     const scalar_type & dt,
                                     Args && ... states) const
  {
    jacobian_type J(3,3);
    this->timeDiscreteJacobian(step, time, dt, J, std::forward<Args>(states)...);
    return J;
  }
};



struct Bdf1Solver
{
  using app_t		= MyApp;
  using nstate_t	= typename app_t::state_type;
  using nveloc_t	= typename app_t::velocity_type;
  using njacobian_t	= typename app_t::jacobian_type;

  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nveloc_t>;
  using jac_t		= ::pressio::containers::Matrix<njacobian_t>;

  using stepper_t = ::pressio::ode::ImplicitStepper<::pressio::ode::ImplicitEnum::Euler,
						    state_t, res_t, jac_t, app_t>;

  using lin_solver_name = ::pressio::solvers::linear::iterative::Bicgstab;
  using lin_solver_t = ::pressio::solvers::iterative::EigenIterative<lin_solver_name, jac_t>;
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<double, lin_solver_t>;

  state_t y_ = {};
  app_t appObj_ = {};
  stepper_t stepperObj_;

  Bdf1Solver()
    : appObj_{}, stepperObj_{y_, appObj_}
  {
    y_.resize(3);
    y_[0] = 1; y_[1]=2; y_[2]=3;
  }

  void run(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(linSolverObj);
    // integrate in time
    const double dt = 0.01;
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt, steps, solverO);
  };
};


struct CustomBdf1Solver
{
  using app_t		= MyApp;
  using nstate_t	= typename app_t::state_type;
  using nresid_t	= typename app_t::residual_type;
  using njacobian_t	= typename app_t::jacobian_type;

  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nresid_t>;
  using jac_t		= ::pressio::containers::Matrix<njacobian_t>;

  using my_custom_order = ::pressio::ode::types::StepperOrder<1>;
  using my_num_states	= ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<::pressio::ode::ImplicitEnum::Arbitrary,
						    state_t, res_t, jac_t, app_t,
						    my_custom_order, my_num_states>;

  using traits = ::pressio::ode::details::traits<stepper_t>;
  using rp_t = typename traits::residual_policy_t;
  static_assert( !std::is_void<rp_t>::value, "");

  using std_r1 = ::pressio::ode::policy::ImplicitResidualStandardPolicyForArbitraryStepper<state_t, app_t, res_t>;
  using std_r2 = typename traits::standard_res_policy_t;
  static_assert( std::is_same<rp_t, std_r1>::value, "");

  using lin_solver_name = ::pressio::solvers::linear::iterative::Bicgstab;
  using lin_solver_t = ::pressio::solvers::iterative::EigenIterative<lin_solver_name, jac_t>;
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<double, lin_solver_t>;

  state_t y_ = {};
  app_t appObj_ = {};
  stepper_t stepperObj_;

  CustomBdf1Solver()
    : appObj_{}, stepperObj_{y_, appObj_}
  {
    y_.resize(3);
    y_[0] = 1; y_[1]=2; y_[2]=3;
  }

  void run(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(linSolverObj);
    // integrate in time
    const double dt = 0.01;
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt, steps, solverO);
  };
};

TEST(ode_implicit, arbitraryStepperRunEuler)
{
  CustomBdf1Solver S1;
  S1.run(2);
  std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

  Bdf1Solver S2;
  S2.run(2);
  std::cout << std::setprecision(14) << *S2.y_.data() << "\n";


  // appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);
  // EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  // EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  // EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);

}







// template<typename state_type, typename system_type, typename residual_type>
// class ResidualPolicy
//   : public ::pressio::ode::policy::ImplicitResidualPolicyBase<
//   ResidualPolicy<state_type, system_type, residual_type>
//   >
// {
// public:
//   static constexpr auto stepper_order = 1;
//   static constexpr auto num_aux_states = 1;

//   void operator()(const state_type & y,
// 		  residual_type & R,
// 		  const std::array<state_type, num_aux_states> & oldYs,
// 		  const system_type & model,
// 		  double t,
// 		  double dt,
// 		  ::pressio::ode::types::step_t step) const
//   {
//     model.velocity(*y.data(), t, *R.data());
//     R.data()->setConstant(1);
//   }

//   residual_type operator()(const state_type & y,
//   			   const std::array<state_type, num_aux_states> & oldYs,
//   			   const system_type & model,
//   			   double t,
//   			   double dt,
// 			   ::pressio::ode::types::step_t step) const{
//     // here I would need to compute the time discrete residual
//     residual_type R(3);
//     R.data()->setConstant(1);
//     return R;
//   }
// };//end class


// template<typename state_type, typename system_type, typename jacobian_type>
// class JacobianPolicy
//   : public ::pressio::ode::policy::JacobianPolicyBase<
//   JacobianPolicy<state_type, system_type, jacobian_type>
//   >
// {
// public:
//   static constexpr auto stepper_order = 1;
//   static constexpr auto num_aux_states = 1;

//   void operator()(const state_type & y,
// 		  jacobian_type & J,
// 		  const system_type & model,
// 		  double t,
// 		  double dt,
// 		  ::pressio::ode::types::step_t step) const{
//     J.resize(3,3);
//     // here I would need to compute the time discrete version
//   }

//   jacobian_type operator()(const state_type & y,
//   			   const system_type & model,
//   			   double t,
//   			   double dt,
// 			   ::pressio::ode::types::step_t step) const{
//     // here I would need to compute the time discrete version
//     jacobian_type J(3, 3);
//     return J;
//   }

// };//end class
