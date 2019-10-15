
#include <gtest/gtest.h>
#include "ROM_LSPG"
#include "../epetra_skeleton.hpp"

// using namespace pressio;
// using matrix_w_t  = containers::MultiVector<Epetra_MultiVector>;
// using decoder_t = rom::LinearDecoder<matrix_w_t>;
// using fom_state_w_t = containers::Vector<Epetra_Vector>;
// using fom_states = rom::FomStatesData<fom_state_w_t, 1, decoder_t>;
// using rom_state_t = containers::Vector<Eigen::VectorXd>;

// struct mytest : pressio::rom::FomStatesData<fom_state_w_t, 1, decoder_t>{
//   using base_t = pressio::rom::FomStatesData<fom_state_w_t, 1, decoder_t>;

//   mytest(const fom_state_w_t & y0Fom, const decoder_t & decObj)
//     : base_t(y0Fom, decObj){
//     fom_states MyStates(y0Fom, decObj);
//     rom_state_t rY(2);
//     rY.putScalar(1.2);
//     this->reconstructCurrentFomState(rY);
//   }

//   void check(){
//     for (auto i=0; i<this->yFom_.localSize(); i++)
//       EXPECT_DOUBLE_EQ(this->yFom_[i], 2.4);
//   }
// };

TEST(lspg, epetra_types)
{
  using fom_t		= pressio::rom::test::EpetraSkeleton;
  using scalar_t	= typename fom_t::scalar_type;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  // define LSPG type
  using lspg_problem_types = pressio::rom::DefaultLSPGTypeGenerator<
    fom_t, pressio::ode::ImplicitEnum::Euler, decoder_t, lspg_state_t>;

  using lspg_stepper_t = typename lspg_problem_types::lspg_stepper_t;

  static_assert(!std::is_void<lspg_stepper_t>::value, "");
}