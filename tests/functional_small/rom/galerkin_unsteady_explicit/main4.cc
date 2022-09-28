
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"

struct MatrixType{};
struct FomStateType{};
struct FomRhsType{};
struct MaskedRhs{};

struct MyFom
{
  using time_type = double;
  using state_type = FomStateType;
  using right_hand_side_type = FomRhsType;

  right_hand_side_type createRightHandSide() const { return right_hand_side_type{}; }
  void rightHandSide(const state_type & u, time_type, right_hand_side_type & r) const{}
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
template<> struct Traits<MaskedRhs>{
  using scalar_type = double;
  static constexpr int rank = 1;
};
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
#include "pressio/rom_galerkin_unsteady.hpp"

struct MyMasker
{
  using operand_type = FomRhsType;
  using result_type = MaskedRhs;

  result_type createApplyMaskResult(const operand_type & /*operand*/) const{
    return result_type{};
  }
  void operator()(const operand_type & operand, result_type & result) const{}
};

struct MyHypRedOperator
{
  using time_type = double;
  using right_hand_side_operand_type = MaskedRhs;

  template<class ResultType>
  void operator()(const right_hand_side_operand_type & operand,
		  time_type,
		  ResultType & result) const{}
};

TEST(rom_galerkin_unsteady, test4)
{
  /* this test is compile-time only and checks that
     for masked unsteady explicit galerkin all the various
     layers can have/use different types and things behave as they should

     - fom uses FomStateType, FomRhsType
     - the masker acts on the fom types but produces "MaskedRhs"
     - the  hypRedOp acts on the masked types
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  MyFom fomSystem;
  MatrixType phi;

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift;
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  MyMasker m1;
  MyHypRedOperator hrOp;
  namespace pg = pressio::rom::galerkin;
  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  auto problem = pg::create_unsteady_explicit_problem(odeScheme, space, fomSystem, m1, hrOp);

  pressio::ode::advance_n_steps(problem, romState, 0.,
				double{}, ::pressio::ode::StepCount(1));

  pressio::log::finalize();
}
