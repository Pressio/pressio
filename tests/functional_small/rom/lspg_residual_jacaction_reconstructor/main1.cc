
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
  using discrete_residual_type = state_type;

  explicit MyFom(fixture_t & a)
    : contigMap_(a.contigMap_){}

  discrete_residual_type createDiscreteTimeResidual() const{
    vec_t a(contigMap_);
    pressio::ops::fill(a, 0);
    return a;
  }

  phi_type createResultOfDiscreteTimeJacobianActionOn(const phi_type & B) const{
    phi_type a(contigMap_, B.getNumVectors());
    pressio::ops::fill(a, 0);
    return a;
  }

  template<class StepCountType>
  void discreteTimeResidualAndJacobianAction(StepCountType,
					     double time,
					     double dt,
					     discrete_residual_type & R,
					     const phi_type & B,
					     std::optional<phi_type*> JA,
					     const state_type & y_np1,
					     const state_type & y_n ) const
  {
    // totally arbitrary operations

    // R = 2 * y_np1 + y_n * 3
    pressio::ops::update(R, 0, y_np1, 2., y_n, 3);

    if (bool(JA)){
      auto & J = *JA.value();
      auto J_view = J.getLocalViewHost(Tpetra::Access::ReadWrite);

      int count = 0;
      for (std::size_t j=0; j< pressio::ops::extent(J_view, 1); ++j){
	for (std::size_t i=0; i< pressio::ops::extent(J_view, 0); ++i){
	  J_view(i,j) = time + count++;
	}
      }
    }

    auto R_view = R.getLocalViewHost(Tpetra::Access::ReadOnly);

    double goldVal = 0.;
    if (time == 2.) { goldVal = 12.; }
    if (time == 4.) { goldVal = 42.; }
    if (time == 6.) { goldVal = 72.; }
    if (time == 8.) { goldVal = 102.; }
    if (time == 10.){ goldVal = 132.; }
    for (int i=0; i<5; ++i){
      ASSERT_DOUBLE_EQ( R_view(i,0), goldVal );
    }

  }
};

TEST_F(fixture_t, lspg_residual_jacaction_reconstructor)
{
  using namespace pressio::rom;

  auto phi = ::pressio::ops::clone(*myMv_);
  for (int i=0; i<numVecs_; ++i){
    auto col = pressio::column(phi, i);
    pressio::ops::fill(col, i);
  }

  vec_t shift(contigMap_);
  pressio::ops::fill(shift, 0);
  using reduced_state_type = Eigen::VectorXd;
  auto space = create_trial_column_subspace<reduced_state_type>(std::move(phi), shift, false);

  const std::string romStatesStr = "./lspg_residual_jacaction_reconstructor/rom_states.txt";
  MyFom app(*this);
  auto o = lspg::create_reconstructor(space);
  o.execute<2>(app, romStatesStr, "main1_");
}
