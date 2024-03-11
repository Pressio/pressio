
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom.hpp"

namespace
{

namespace pode = pressio::ode;
namespace prom = pressio::rom;

constexpr int _m = 8;
constexpr int _n = 3;
constexpr pode::StepScheme odeScheme = pode::StepScheme::BDF1;

using vec_t = Eigen::VectorXd;
using phi_t = Eigen::Matrix<double, -1,-1>;

void write_vec_cout(const std::string & s, const Eigen::VectorXd & v){
  std::cout << s << " ";
  for (int i=0; i<v.size(); ++i){
    std::cout << v[i] << " ";
  }
}

struct MyFom
{
  using independent_variable_type = double;
  using time_type  = double;
  using state_type = vec_t;
  using rhs_type   = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type createState() const{
    state_type s(_m);
    s.setConstant(0.);
    return s;
  }

  rhs_type createRhs() const{
    rhs_type R(_m);
    R.setConstant(0.);
    return R;
  }

  jacobian_type createJacobian() const{
    jacobian_type JJ(_m,_m);
    return JJ;
  };

  template<class OperandType>
  OperandType createResultOfJacobianActionOn(const OperandType & B) const{
    OperandType A(_m, B.cols());
    A.setConstant(0.);
    return A;
  }

  void rhsAndJacobian(const state_type & yIn,
		      independent_variable_type /*unused*/,
		      rhs_type & R,
#ifdef PRESSIO_ENABLE_CXX17
		      std::optional<jacobian_type*> Jin) const
#else
                      jacobian_type* Jin) const
#endif
  {
    R.setConstant(1.);
  }

  void rhs(const state_type & yIn,
	   independent_variable_type /*unused*/,
	   rhs_type & R) const
  {
    R.setConstant(1.);
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     time_type time,
                     OperandType & A) const
  {}
};

struct MyFakeNonLinSolverForFOM{
  int count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    count_++;

    using r_t = typename SystemType::residual_type;
    r_t goldR3(_m); goldR3.setConstant(-2.);
    r_t goldR4(_m); goldR4.setConstant(-1.);

    auto R = system.createResidual();
    auto J = system.createJacobian();

    // mimic solve step 1
    system.residualAndJacobian(state, R, &J);
    write_vec_cout("R = ", R); std::cout << "\n";
    if      (count_ == 1){ ASSERT_TRUE( R.isApprox(goldR3) ); }
    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

    // mimic solve step 2
    system.residualAndJacobian(state, R, &J);
    write_vec_cout("R = ", R); std::cout << "\n";
    if      (count_ == 1){ ASSERT_TRUE( R.isApprox(goldR4) ); }
    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
  }
};

struct MyFakeNonLinSolverForROM{
  int count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    count_++;

    using r_t = typename SystemType::residual_type;
    r_t goldR1(_m); goldR1.setConstant(-2.);
    r_t goldR2(_m); goldR2.setConstant(4.);
    r_t goldR5(_m); goldR5.setConstant(-2.);
    r_t goldR6(_m); goldR6.setConstant(4.);

    auto R = system.createResidual();
    auto J = system.createJacobian();

    // mimic solve step 1
    system.residualAndJacobian(state, R, &J);
    write_vec_cout("R = ", R); std::cout << "\n";
    if      (count_ == 1){ ASSERT_TRUE( R.isApprox(goldR1) ); }
    else if (count_ == 2){ ASSERT_TRUE( R.isApprox(goldR5) ); }
    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

    // mimic solve step 2
    system.residualAndJacobian(state, R, &J);
    write_vec_cout("R = ", R); std::cout << "\n";
    if      (count_ == 1){ ASSERT_TRUE( R.isApprox(goldR2) ); }
    else if (count_ == 2){ ASSERT_TRUE( R.isApprox(goldR6) ); }
    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
  }
};

struct Subdomain
{
  using fom_state_t   = typename MyFom::state_type;
  using rom_state_t   = Eigen::VectorXd;
  using trial_subspace_t =
    decltype(prom::create_trial_column_subspace<rom_state_t>(std::declval<const phi_t &>(),
							     std::declval<fom_state_t&&>(),
							     false));
  using mixed_problem_t =
    decltype(prom::lspg::experimental::create_unsteady_problem_mixed_fom(odeScheme,
									 std::declval<trial_subspace_t const &>(),
									 std::declval<MyFom const &>()));

  Subdomain(const vec_t & shift, const phi_t & phiIn)
    : fomObj_(), fomState_(_m), romState_(_n),
      romSubspace_(prom::create_trial_column_subspace<rom_state_t>(phiIn, shift, false)),
      mixedProblem_(prom::lspg::experimental::create_unsteady_problem_mixed_fom(odeScheme,
										romSubspace_,
										fomObj_)),
      fomNonLinSolver_(),
      romNonLinSolver_()
  {
    fomState_.setConstant(10.);
    romState_.setConstant(2.);
  }

  void doStep(const double startTime, const int step, const double dtIn)
  {
    const pode::StepStartAt<double> startAt(startTime);
    const pode::StepCount stepNumber(step);
    const pode::StepSize<double> dt(dtIn);
    std::cout << "\n***** DOING STEP " << step << " *****\n";
    std::cout << "-----------------------------\n";

    write_vec_cout("before step, fomState = ", fomState_); std::cout << '\n';
    write_vec_cout("before step, romState = ", romState_); std::cout << '\n';

    if (step == 1){
      ASSERT_DOUBLE_EQ(startTime, 0.);
      mixedProblem_(prom::RomStep, romState_, startAt, stepNumber, dt, romNonLinSolver_);
    }

    else if (step == 2){
      ASSERT_DOUBLE_EQ(startTime, 0.+dtIn);
      mixedProblem_(prom::TransitionToFomAndDoStep, fomState_, romState_, startAt,
		    stepNumber, dt, fomNonLinSolver_);
    }
    else if (step == 3){
      ASSERT_DOUBLE_EQ(startTime, 2*dtIn);
      mixedProblem_(prom::TransitionToRomAndDoStep, fomState_, romState_, startAt,
		    stepNumber, dt, romNonLinSolver_);
    }

    write_vec_cout("after step, fomState = ", fomState_); std::cout << '\n';
    write_vec_cout("after step, romState = ", romState_); std::cout << '\n';
  }

  int fomSolverCount() const { return fomNonLinSolver_.count_; }
  int romSolverCount() const { return romNonLinSolver_.count_; }
  auto fomState() const{ return fomState_; }
  auto romState() const{ return romState_; }

private:
  fom_state_t getFomState() const{ return fomState_; }

  MyFom fomObj_;
  fom_state_t fomState_;
  rom_state_t romState_;
  trial_subspace_t romSubspace_;
  mixed_problem_t mixedProblem_;
  MyFakeNonLinSolverForFOM fomNonLinSolver_;
  MyFakeNonLinSolverForROM romNonLinSolver_;
};
}

TEST(rom_lspg_unsteady, test_bdf1_rom_fom_rom)
{
  /*
  stepID     1          2          3

       o -------- o -------- o -------- o
           rom        fom         rom
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  vec_t shift(_m);
  shift.setConstant(0.);

  phi_t phi(_m, _n);
  for (int j=0; j<_n; ++j){
    phi.col(j).setConstant( (double) (j+1));
  }
  std::cout << phi << "\n";

  Subdomain subdom(shift, phi);

  const int numFomSteps = 1;
  const int numRomSteps = 2;
  const int numSteps = numFomSteps + numRomSteps;
  const double dt = 2.;
  double time = 0.;
  for (int i=1; i<=numSteps; ++i){
    subdom.doStep(time, i, dt);
    time += dt;
  }

  // verify things
  ASSERT_TRUE(subdom.fomSolverCount() == numFomSteps);
  ASSERT_TRUE(subdom.romSolverCount() == numRomSteps);

  const auto & computedFomState = subdom.fomState();
  const auto & computedRomState = subdom.romState();

  Eigen::VectorXd goldFomState(_m);
  pressio::ops::fill(goldFomState, 26.);
  ASSERT_TRUE( goldFomState.isApprox(computedFomState) );
  Eigen::VectorXd goldRomState(_n);
  goldRomState[0] = 210.; goldRomState[1] = 418.; goldRomState[2] = 626.;
  ASSERT_TRUE( goldRomState.isApprox(computedRomState) ) << computedRomState;

  pressio::log::finalize();
}
