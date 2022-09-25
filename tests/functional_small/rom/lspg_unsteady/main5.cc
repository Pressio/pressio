
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"

struct MatrixType{};
struct FomStateType{};
struct FomRhsType{};
struct FomJacobianActionResultType{};

struct MyFom
{
  using time_type = double;
  using state_type = FomStateType;
  using right_hand_side_type = FomRhsType;

  right_hand_side_type createRightHandSide() const{ return right_hand_side_type{}; }
  FomJacobianActionResultType createApplyJacobianResult(const MatrixType & B) const{
    return {};
  }
  void rightHandSide(const state_type & u,
		     time_type timeIn,
		     right_hand_side_type & r) const{}
  void applyJacobian(const state_type & state,
                     const MatrixType & B,
                     time_type time,
		     FomJacobianActionResultType) const{}
};

namespace pressio{
template<> struct Traits<MatrixType>{
  using scalar_type = double;
  static constexpr int rank = 2;
};
template<> struct Traits<FomStateType>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
template<> struct Traits<FomRhsType>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
template<> struct Traits<FomJacobianActionResultType>{
  using scalar_type = double;
  static constexpr int rank = 2;
};
}

#include "pressio/ops.hpp"

namespace pressio { namespace ops{

//
// needed by trial subspace
//
std::size_t extent(const MatrixType &, int ){ return std::size_t{}; }
MatrixType clone(const MatrixType &){ return MatrixType{}; }

template<class AlphaType, class BetaType>
void product(::pressio::nontranspose /*mode*/,
             const AlphaType & alpha,
             const MatrixType & basis,
             const Eigen::VectorXd & operand,
             const BetaType & beta,
             FomStateType & fullState){}

FomStateType clone(const FomStateType &){ return FomStateType{}; }
template<class ScalarType> void fill(FomStateType & fullState, ScalarType value){}
template<class ScalarType>
void update(FomStateType &, const ScalarType &,
            const FomStateType &, const ScalarType &){}

//
// needed by fom states manager
void deep_copy(FomStateType &, const FomStateType &){}
void set_zero(FomStateType &){}

//
// needed by ode rhs manager
void set_zero(FomRhsType &){}

//
// needed by lspg impl
template<class AlphaType, class BetaType>
void update(FomJacobianActionResultType &, const AlphaType &,
            const MatrixType &, const BetaType&){}

template<class at, class bt, class ct>
void update(FomRhsType &, const at &,
            const FomStateType &, const bt&,
	    const FomStateType &, const ct&){}

template<class at, class bt, class ct, class dt>
void update(FomRhsType &, const at &,
            const FomStateType &, const bt&,
	    const FomStateType &, const ct&,
	    const FomStateType &, const dt&){}

}} // end namespace pressio::ops

#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"

struct FakeNonLinSolver
{
  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    // detect things to make sure system has the right API
    static_assert(std::is_same<StateType,
		  decltype(system.createState())>::value, "");
    static_assert(!std::is_void<decltype(system.createResidual())>::value, "");
    static_assert(!std::is_void<decltype(system.createJacobian())>::value, "");

    using r_t = typename SystemType::residual_type;
    using j_t = typename SystemType::jacobian_type;
    static_assert(std::is_void<
		  decltype(system.residualAndJacobian(state,
		  std::declval<r_t&>(),
		  std::declval<j_t&>(),
		  true))
		  >::value, "");
  }
};

TEST(rom_lspg_unsteady, test6)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  MyFom fomSystem;
  MatrixType phi;

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift;
  auto space = pressio::rom::create_trial_column_subspace<
    reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();

  auto problem = pressio::rom::lspg::create_unsteady_problem
    (pressio::ode::StepScheme::BDF1, space, fomSystem);

  FakeNonLinSolver nonLinSolver;
  pressio::ode::advance_n_steps(problem, romState, 0., 1.,
				::pressio::ode::StepCount(1),
				nonLinSolver);

  pressio::log::finalize();
}
