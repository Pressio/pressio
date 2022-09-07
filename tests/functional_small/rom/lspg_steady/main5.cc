#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"

struct MatrixType{};
struct FomStateType{};
struct FomResidualType{};
struct FomJacobianActionResultType{};
struct MaskedResidual{};
struct MaskedJacobianAction{};

struct MyFom
{
  using state_type = FomStateType;
  using residual_type = FomResidualType;

  residual_type createResidual() const { return residual_type{}; }

  FomJacobianActionResultType createApplyJacobianResult(const MatrixType & B) const{
    return FomJacobianActionResultType{};
  }

  void residual(const state_type & u,
		residual_type & r) const{}

  void applyJacobian(const state_type & state,
                     const MatrixType & B,
		     FomJacobianActionResultType & A) const{}
};

namespace pressio{
template<> struct Traits<MatrixType>{ using scalar_type = double; };
template<> struct Traits<FomStateType>{ using scalar_type = double; };
template<> struct Traits<FomResidualType>{ using scalar_type = double; };
template<> struct Traits<FomJacobianActionResultType>{ using scalar_type = double; };
template<> struct Traits<MaskedResidual>{ using scalar_type = double; };
template<> struct Traits<MaskedJacobianAction>{ using scalar_type = double; };
}

namespace pressio { namespace ops{

std::size_t extent(const MatrixType &, int ){ return std::size_t{}; }

FomStateType clone(const FomStateType &){ return FomStateType{}; }
MatrixType clone(const MatrixType &){ return MatrixType{}; }

template<class ScalarType>
void fill(FomStateType & fullState, ScalarType value){}

template<class ScalarType>
void update(FomStateType & y, const ScalarType & alpha,
            const FomStateType & x, const ScalarType & beta){}

template<class AlphaType, class BetaType>
void product(::pressio::nontranspose /*mode*/,
             const AlphaType & alpha,
             const MatrixType & basis,
             const Eigen::VectorXd & operand,
             const BetaType & beta,
             FomStateType & fullState){}

}} // end namespace pressio::ops

#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_steady.hpp"

struct MyMaskerResidual
{
  using operand_type = FomResidualType;
  using result_type = MaskedResidual;

  result_type createApplyMaskResult(const operand_type & /*operand*/) const{
    return result_type{};
  }

  void operator()(const operand_type & operand, result_type & result) const{}
};

struct MyMaskerJacAction
{
  using operand_type = FomJacobianActionResultType;
  using result_type = MaskedJacobianAction;

  result_type createApplyMaskResult(const operand_type & /*operand*/) const{
    return result_type{};
  }

  void operator()(const operand_type & operand, result_type & result) const{}
};

struct MyHypRedOperator
{
  using residual_operand_type = MaskedResidual;
  using jacobian_action_operand_type = MaskedJacobianAction;

  template<class ResultType>
  void operator()(const residual_operand_type & operand,
		  ResultType & result) const{}

  template<class ResultType>
  void operator()(const jacobian_action_operand_type & operand,
		  ResultType & result) const{}
};

struct FakeNonLinSolver
{
  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    auto R = system.createResidual();
    auto J = system.createJacobian();
    system.residualAndJacobian(state, R, J, true);
  }
};

TEST(rom_lspg_steady, test5)
{
  /* this test is compile-time only and checks that
     for masked steady lspg all the various
     layers can have/use different types and things behave as they should

     - fom uses FomStateType, FomResidualType, FomJacobianActionResultType
     - the masker acts on the fom types but produces "MaskedResidual and MaskedJacobianAction"
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  MyFom fomSystem;
  MatrixType phi;

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift;
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  MyHypRedOperator proj;
  MyMaskerResidual  m1;
  MyMaskerJacAction m2;
  auto problem = pressio::rom::lspg::create_steady_problem(space, fomSystem,
							   m1, m2);

  FakeNonLinSolver nonLinSolver;
  nonLinSolver.solve(problem, romState);

  pressio::log::finalize();
}
