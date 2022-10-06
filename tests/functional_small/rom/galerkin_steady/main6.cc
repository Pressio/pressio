
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"

struct MatrixType{};
struct FomStateType{};
struct FomResidualType{};
struct FomJacobianActionResultType{};

struct MyFom
{
  using state_type = FomStateType;
  using residual_type = FomResidualType;

  residual_type createResidual() const { return residual_type{}; }
  FomJacobianActionResultType createResultOfJacobianActionOn(const MatrixType & B) const{
    return FomJacobianActionResultType{};
  }
  void residual(const state_type & u, residual_type & r) const{}
  void residualAndJacobianAction(const state_type & state,
				 residual_type & r, const MatrixType & B,
				 FomJacobianActionResultType & A,
				 bool /*computeJac*/) const{}
};

namespace pressio{
template<> struct Traits<MatrixType>{
  static constexpr int rank = 2;
  using scalar_type = double;
};

template<> struct Traits<FomStateType>{
  static constexpr int rank = 1;
  using scalar_type = double;
};

template<> struct Traits<FomResidualType>{
  static constexpr int rank = 1;
  using scalar_type = double;
};

template<> struct Traits<FomJacobianActionResultType>{
  static constexpr int rank = 2;
  using scalar_type = double;
};
}

namespace pressio { namespace ops{

template<class AlphaType, class BetaType>
void product(::pressio::transpose /*mode*/,
             const AlphaType & /*alpha*/,
             const MatrixType & /*basis*/,
             const FomResidualType & /*operand*/,
             const BetaType & /*beta*/,
	     Eigen::VectorXd &){}

template<class AlphaType, class BetaType>
void product(::pressio::transpose /*mode*/,
	     ::pressio::nontranspose /*mode*/,
             const AlphaType & /*alpha*/,
             const MatrixType & /*basis*/,
             const FomJacobianActionResultType & /*operand*/,
             const BetaType & /*beta*/,
             Eigen::MatrixXd & /*reducedJac*/){}

}} // end namespace pressio::ops

// this include has to be here because the product
// specializations above have to be visible
#include "pressio/rom_galerkin_steady.hpp"

class MyFakeTrialSubspace
{
  FomStateType fs_;
  MatrixType M_;
  // must be here to satisfy the concept
  MyFakeTrialSubspace & operator=(const MyFakeTrialSubspace &) = delete;

public:
  using reduced_state_type = Eigen::VectorXd;
  using full_state_type    = FomStateType;
  using basis_matrix_type  = MatrixType;

  reduced_state_type createReducedState() const{ return {};}
  full_state_type    createFullState() const{ return fs_;}
  void               mapFromReducedState(const reduced_state_type &, full_state_type &) const{}
  full_state_type    createFullStateFromReducedState(const reduced_state_type &) const{ return fs_;}
  const basis_matrix_type & basisOfTranslatedSpace() const{ return M_;}
  const full_state_type   & translationVector() const{ return fs_;}
  const basis_matrix_type & basis() const{ return M_;}
  std::size_t dimension() const{ return {}; }
  bool isColumnSpace() const{ return true; }
  bool isRowSpace() const{ return false; }
};

struct FakeNonLinSolver{
  template<class SystemType, class StateType>
  void solve(const SystemType & system,
	     StateType & state)
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

TEST(rom_galerkin_steady, test6)
{
  /* this test is compile-time only and is meant to
     check the scenario of usinga custom trial subspace
     that should follow the related concept and
     we can verify the actual operations and things needed
     are those shown in the concept
  */

  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  MyFom fomSystem ;
  MyFakeTrialSubspace space;
  namespace pg = pressio::rom::galerkin;
  auto problem = pg::create_steady_problem(space, fomSystem);
  FakeNonLinSolver nonLinSolver;
  MyFakeTrialSubspace::reduced_state_type romState(5);
  nonLinSolver.solve(problem, romState);

  pressio::log::finalize();
}
