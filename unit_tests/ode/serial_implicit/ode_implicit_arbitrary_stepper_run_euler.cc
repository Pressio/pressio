
#include <gtest/gtest.h>
#include "ODE_ALL"


template <typename state_t>
struct Observer{
  using matrix_t = Eigen::MatrixXd;
  using step_t = pressio::ode::types::step_t;

  size_t state_size_ {};
  matrix_t A_;
  size_t count_ {};

  Observer(int N, int state_size)
    : A_(state_size, N+1){} //+1 to store also init cond

  void operator()(step_t step, double t, const state_t & y)
  {
    this->storeInColumn(y, count_);
    count_++;
  }

  void storeInColumn(const state_t & y, int j){
    for (auto i=0; i<y.size(); i++)
      A_(i,j) = y(i);
  }
  size_t getCount() const{ return count_;}
  void printAll() const{ std::cout << A_ << std::endl; }

  const matrix_t & viewSnaps() const{
    return A_;
  }
};


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

public:
  void velocity(const state_type & yIn,
		const scalar_type & t,
		velocity_type & f) const{
    f = -10. * yIn;
  };

  velocity_type velocity(const state_type & yIn,
			 const scalar_type & t) const{
    velocity_type f(3);
    this->velocity(yIn, t, f);
    return f;
  };

  void jacobian(const state_type & yIn,
  		const scalar_type & t,
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
  			 const scalar_type & t) const{
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

  residual_type createTimeDiscreteResidualObject(const state_type & state) const
  {
    residual_type R(3);
    R.setConstant(0);
    return R;
  }

  jacobian_type createTimeDiscreteJacobianObject(const state_type & state) const
  {
    jacobian_type J(3,3);
    typedef Eigen::Triplet<scalar_type> Tr;
    std::vector<Tr> tripletList;
    tripletList.push_back( Tr( 0, 0, 0.) );
    tripletList.push_back( Tr( 1, 1, 0.) );
    tripletList.push_back( Tr( 2, 2, 0.) );
    J.setFromTriplets(tripletList.begin(), tripletList.end());
    return J;
  }

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

};



struct Bdf1Solver
{
  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
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
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, sc_t>;

  app_t appObj_ = {};
  state_t y_ = {};
  stepper_t stepperObj_;
  const sc_t dt_ = 0.01;

  Bdf1Solver(const state_t & yIn)
    : appObj_{}, y_{yIn}, stepperObj_{y_, appObj_}
  {}

  void integrateForNSteps(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt_, steps, solverO);
  };
};


struct CustomBdf1Solver
{
  using app_t		= MyApp;
  using sc_t		= typename app_t::scalar_type;
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

  using std_r1 = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicyForArbitraryStepper<state_t, app_t, res_t>;
  using std_r2 = typename traits::standard_res_policy_t;
  static_assert( std::is_same<rp_t, std_r1>::value, "");

  using lin_solver_name = ::pressio::solvers::linear::iterative::Bicgstab;
  using lin_solver_t = ::pressio::solvers::iterative::EigenIterative<lin_solver_name, jac_t>;
  using nonlin_solver_t = ::pressio::solvers::NewtonRaphson<stepper_t, lin_solver_t, sc_t>;

  app_t appObj_		= {};
  state_t y_		= {};
  stepper_t stepperObj_;
  const sc_t dt_	= 0.01;

  CustomBdf1Solver(const state_t & yIn)
    : appObj_{}, y_{yIn}, stepperObj_{y_, appObj_}
  {}

  void integrateForNSteps(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, dt_, steps, solverO);
  };

  void integrateForNStepsWithStepSizeManagerLambda(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_;
				};
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, steps, solverO, dtSetterLambda);
  };

  void integrateForNStepsWithStepSizeManagerLambdaWrongDt(int steps)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_*2.;
				};
    ::pressio::ode::integrateNSteps(stepperObj_, y_, 0.0, steps, solverO, dtSetterLambda);
  };

  void integrateToTimeWithStepSizeManagerLambda(double finalTime)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_;
				};
    ::pressio::ode::integrateToTargetTime(stepperObj_, y_, 0.0, finalTime, solverO, dtSetterLambda);
  };


  template <typename observer_t>
  void integrateToTimeWithStepSizeManagerLambdaAndCollector(double finalTime, observer_t & observer)
  {
    lin_solver_t linSolverObj;
    nonlin_solver_t solverO(stepperObj_, y_, linSolverObj);
    using step_t = ::pressio::ode::types::step_t;
    const auto dtSetterLambda = [=](const step_t & step, const sc_t & time, sc_t & dt){
				  std::cout << " SETTING DT " << std::endl;
				  dt = dt_;
				};
    ::pressio::ode::integrateToTargetTime(stepperObj_, y_, 0.0,
					  finalTime, solverO,
					  dtSetterLambda, observer);
  };

};



TEST(ode_implicit, arbitraryStepperRunEulerConstDt)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    S1.integrateForNSteps(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.integrateForNSteps(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    EXPECT_DOUBLE_EQ( S1.y_[0], S2.y_[0]);
    EXPECT_DOUBLE_EQ( S1.y_[1], S2.y_[1]);
    EXPECT_DOUBLE_EQ( S1.y_[2], S2.y_[2]);
  }
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetter)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    S1.integrateForNStepsWithStepSizeManagerLambda(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.integrateForNSteps(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    EXPECT_DOUBLE_EQ( S1.y_[0], S2.y_[0]);
    EXPECT_DOUBLE_EQ( S1.y_[1], S2.y_[1]);
    EXPECT_DOUBLE_EQ( S1.y_[2], S2.y_[2]);
  }
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetterWithWrongDt)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  for (int N = 1; N < 10; N++){
    CustomBdf1Solver S1(y0);
    // use here a dt that we know wont work because it is wrong
    S1.integrateForNStepsWithStepSizeManagerLambdaWrongDt(N);
    std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

    Bdf1Solver S2(y0);
    S2.integrateForNSteps(N);
    std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

    ASSERT_TRUE( S1.y_[0] != S2.y_[0]);
    ASSERT_TRUE( S1.y_[1] != S2.y_[1]);
    ASSERT_TRUE( S1.y_[2] != S2.y_[2]);
  }
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetterIntegrateToTimeTrivial)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  // here we set final time to zero so this should not do anything,
  // the state should remain the same as init condition
  constexpr double finalTime = 0.0;

  CustomBdf1Solver S1(y0);
  S1.integrateToTimeWithStepSizeManagerLambda(finalTime);
  std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

  EXPECT_DOUBLE_EQ( S1.y_[0] , one );
  EXPECT_DOUBLE_EQ( S1.y_[1] , two );
  EXPECT_DOUBLE_EQ( S1.y_[2] , three);
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetterIntegrateToTimeNonTrivial)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  // here we set final time to something non trivial and such that
  // it is divisble by the default dt so that we get an integer number of steps
  constexpr double finalTime = 0.1;

  CustomBdf1Solver S1(y0);
  S1.integrateToTimeWithStepSizeManagerLambda(finalTime);
  std::cout << std::setprecision(14) << *S1.y_.data() << "\n";

  Bdf1Solver S2(y0);
  // compute num of steps so that we make sure the integrate
  // to target time matches the integrate n steps
  const auto dt = S2.dt_;
  const int nSteps = finalTime/dt;
  S2.integrateForNSteps(nSteps);
  std::cout << std::setprecision(14) << *S2.y_.data() << "\n";

  EXPECT_DOUBLE_EQ( S1.y_[0], S2.y_[0]);
  EXPECT_DOUBLE_EQ( S1.y_[1], S2.y_[1]);
  EXPECT_DOUBLE_EQ( S1.y_[2], S2.y_[2]);
}


TEST(ode_implicit, arbitraryStepperRunEulerDtSetterIntegrateToTimeNonTrivialWithCollector)
{
  constexpr double one = ::pressio::utils::constants::one<double>();
  constexpr double two = ::pressio::utils::constants::two<double>();
  constexpr double three = ::pressio::utils::constants::three<double>();
  ::pressio::containers::Vector<Eigen::VectorXd> y0(3);
  *y0.data() << one,two,three;

  // here we set final time to something non trivial and such that
  // it is divisble by the default dt so that we get an integer number of steps
  constexpr double finalTime = 0.1;

  // define observer (templated on the native state type so that pressio
  // passes us an object of the native state)
  using state_t = typename MyApp::state_type;
  using observer_t = Observer<state_t>;

  // create test class object
  CustomBdf1Solver S1(y0);

  const auto dt = S1.dt_;
  const int nSteps = finalTime/dt;
  // create observer
  observer_t observerObj(nSteps, y0.extent(0));

  S1.integrateToTimeWithStepSizeManagerLambdaAndCollector(finalTime, observerObj);
  std::cout << std::setprecision(14) << *S1.y_.data() << "\n";
  observerObj.printAll();

  // check that observer has right data
  Eigen::MatrixXd trueS(y0.extent(0), nSteps+1);
  trueS(0,0)=1.00000000000000e+00;
  trueS(1,0)=2.00000000000000e+00;
  trueS(2,0)=3.00000000000000e+00;

  trueS(0,2)=8.26446280991735e-01;
  trueS(1,2)=1.65289256198347e+00;
  trueS(2,2)=2.47933884297521e+00;

  trueS(0,6)=5.64473930053777e-01;
  trueS(1,6)=1.12894786010755e+00;
  trueS(2,6)=1.69342179016133e+00;

  trueS(0,10)=3.85543289429532e-01;
  trueS(1,10)=7.71086578859064e-01;
  trueS(2,10)=1.15662986828860e+00;

  const auto & storedData = observerObj.viewSnaps();
  std::vector<int> colInd = {0,2,6,10};
  for (const auto & j : colInd){
    EXPECT_NEAR( storedData(0,j), trueS(0,j), 1e-13);
    EXPECT_NEAR( storedData(1,j), trueS(1,j), 1e-13);
    EXPECT_NEAR( storedData(2,j), trueS(2,j), 1e-13);
  }
}
