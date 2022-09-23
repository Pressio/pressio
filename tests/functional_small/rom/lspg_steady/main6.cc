#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/ops.hpp"
#include "pressio/rom_lspg_steady.hpp"

struct MatrixType{};
struct FomStateType{};
struct FomResidualType{};
struct FomJacobianActionResultType{};

struct MyFom
{
  using state_type = FomStateType;
  using residual_type = FomResidualType;

  residual_type createResidual() const { return residual_type{}; }
  FomJacobianActionResultType createApplyJacobianResult(const MatrixType & B) const{
    return FomJacobianActionResultType{};
  }
  void residual(const state_type & u, residual_type & r) const{}
  void applyJacobian(const state_type & state,
                     const MatrixType & B,
		     FomJacobianActionResultType & A) const{}
};

namespace pressio{
template<> struct Traits<MatrixType>{ using scalar_type = double; };
template<> struct Traits<FomStateType>{ using scalar_type = double; };
template<> struct Traits<FomResidualType>{ using scalar_type = double; };
template<> struct Traits<FomJacobianActionResultType>{ using scalar_type = double; };
}

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

TEST(rom_lspg_steady, test6)
{
  /* this test is compile-time only and checks that
     for steady lspg we can provide a custom trial subspace
     that satisfies concept
  */


  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  MyFom fomSystem;
  MyFakeTrialSubspace space;
  auto problem = pressio::rom::lspg::create_steady_problem(space, fomSystem);
  FakeNonLinSolver nonLinSolver;
  MyFakeTrialSubspace::reduced_state_type romState(5);
  nonLinSolver.solve(problem, romState);

  pressio::log::finalize();
}
