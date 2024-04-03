
#include <gtest/gtest.h>
#include "pressio/type_traits.hpp"
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"
#include "../fixtures/tpetra_only_fixtures.hpp"

using fixture_t = tpetraMultiVectorGlobSize15Fixture;

class MyFom
{
  using vec_t = typename fixture_t::vec_t;
  using phi_type = typename fixture_t::mvec_t;
  using map_t = typename fixture_t::map_t;
  Teuchos::RCP<const map_t> contigMap_;

public:
  using time_type = double;
  using state_type = vec_t;
  using rhs_type  = state_type;

  explicit MyFom(fixture_t & a)
    : contigMap_(a.contigMap_){}

  rhs_type createRhs() const{
    vec_t a(contigMap_);
    pressio::ops::fill(a, 0);
    return a;
  }

  phi_type createResultOfJacobianActionOn(const phi_type & B) const{
    phi_type a(contigMap_, B.getNumVectors());
    pressio::ops::fill(a, 0);
    return a;
  }

  void rhs(const state_type & y,
	   double time,
	   rhs_type & f) const
  {
    // totally arbitrary operations
    pressio::ops::fill(f, 1. + time);
  }

  void applyJacobian(const state_type & state,
                     const phi_type & B,
                     time_type time,
		     phi_type & A) const
  {
    auto J_view = A.getLocalViewHost(Tpetra::Access::ReadWrite);
    int count = 0;
    for (std::size_t j=0; j< pressio::ops::extent(J_view, 1); ++j){
      for (std::size_t i=0; i< pressio::ops::extent(J_view, 0); ++i){
	J_view(i,j) = time + count++;
      }
    }
  }
};

TEST_F(fixture_t, lspg_residual_jacaction_reconstructor_bdf1)
{
  using namespace pressio;

  auto phi = ops::clone(*myMv_);
  for (int i=0; i<numVecs_; ++i){
    auto col = pressio::column(phi, i);
    ops::fill(col, i);
  }

  vec_t shift(contigMap_);
  ops::fill(shift, 0);
  using reduced_state_type = Eigen::VectorXd;
  auto space = rom::create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, false);

  const std::string romStatesStr = "./lspg_residual_jacaction_reconstructor/rom_states.txt";
  MyFom app(*this);
  auto o = rom::lspg::create_reconstructor(space);
  o.execute(app, romStatesStr, ode::StepScheme::BDF1, "main2_");
}
