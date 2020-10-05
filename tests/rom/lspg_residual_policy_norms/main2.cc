
#include "pressio_rom.hpp"

struct MyFakeApp
{
  int N_ = {};
public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N) : N_(N){}

  residual_type createResidual() const{ return residual_type(N_);}

  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    dense_matrix_type A;
    return A;
  }

  void residual(const state_type & state,
		residual_type & f) const
  {
    for (int i=0; i<N_; ++i) f(i) = 1;
  }
  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     dense_matrix_type & A) const
  {}

  void applyPreconditioner(const state_type & yState,
			   dense_matrix_type & C) const
  {
  }
  void applyPreconditioner(const state_type & yState,
			   residual_type & operand) const
  {
    for (int i=0; i<N_; ++i) operand(i) += 0.5;
  }

  // residual_type createApplyMaskResult(const residual_type & operand) const{
  //   return residual_type(N_);
  // }
  // dense_matrix_type createApplyMaskResult(const dense_matrix_type & operand) const
  // {
  //   dense_matrix_type A;
  //   return A;
  // }

  // void applyMask(const residual_type & operand,
		//  residual_type & result) const
  // {
  //   for (int i=0; i<N_; ++i) result(i) = operand(i) + 1.5;
  // }

  // void applyMask(const dense_matrix_type &,
		//  dense_matrix_type &) const{}

};


struct MyFakeSolver
{
  MyFakeSolver(){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state){}
};


constexpr int fomSize = 10;
constexpr int romSize = 4;

using scalar_t		= double;
using native_state_t	= Eigen::VectorXd;
using native_dmat_t	= Eigen::MatrixXd;
using fom_state_t	= pressio::containers::Vector<native_state_t>;
using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

using lspg_state_t    = pressio::containers::Vector<Eigen::VectorXd>;
using lspg_residual_t = fom_state_t;


void testSteadyResidualPolicy(fom_state_t & yRef,
			      decoder_t & decoderObj)
{
  std::string checkStr {"PASSED"};

  MyFakeApp fomObj(fomSize);

  // ------------------------------------
  // construct what we need to construct
  // the residual policy object
  // ------------------------------------
  using fom_recon = pressio::rom::FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;
  fom_recon fomReconstructor(yRef, decoderObj);

  using fom_states_manager_t
    = pressio::rom::ManagerFomStatesStatic<fom_state_t, 1, fom_recon, void>;
  fom_states_manager_t fomStatesMngr(fomReconstructor, yRef);

  using lspg_residual_policy
    = ::pressio::rom::lspg::impl::steady::ResidualPolicy<
      lspg_residual_t, fom_states_manager_t, void>;

  using lspg_precond_residual_policy =
    ::pressio::rom::decorator::PreconditionedResidualPolicy<lspg_residual_policy>;

  const lspg_residual_policy rPol(fomStatesMngr);
  const lspg_precond_residual_policy rPolPrecond(fomStatesMngr);

  // ------------------------------------
  // ***** here we do the test *****
  // ------------------------------------
  lspg_state_t romState(romSize);
  lspg_residual_t residual(fomSize);

  // set romY_n = {1,2,3,4}
  romState[0] = 1.;
  romState[1] = 2.;
  romState[2] = 3.;
  romState[3] = 4.;

  /* the compute call below should first compute the residual for BDF1:
  	R = f( phi*romY_n)
     and then call:
  	applyPreconditioned(R) on the fomObj

     where we set f= to be all 1, dt = .5, phi is all 1.
     - apply precon takes R and mimics a modification by adding 0.5 to each entry

     here we don't care if this is meaningful, we just want to test things.
     so we should have:
     fomY_n = [1 1 1 1|   1           10
  	      |1 1 1 1|   2           10
  	      |1 1 1 1|   3     =     10
  	      |1 1 1 1|   4           10
  	      |...	              ...
  	      |1 1 1 1]		      10

     so before preconditioning, we have
  	R  = [1 1 ... 1]

     after precond, we have:
  	R  = [1.5 1.5 ... 1.5]
  */

  scalar_t rNorm1 = {};
  scalar_t rNorm2 = {};

  rPol.compute(romState, residual, fomObj);
  rNorm1 = residual.data()->norm();
  rPolPrecond.compute(romState, residual, fomObj);
  rNorm2 = residual.data()->norm();

  const auto trueNorm1 = std::sqrt( (1*1)*(double)fomSize );
  const auto trueNorm2 = std::sqrt( (1.5*1.5)*(double)fomSize );

  const auto err1 = std::abs( rNorm1 - trueNorm1 );
  const auto err2 = std::abs( rNorm2 - trueNorm2 );

  if ( err1 > 1e-12 or err2 > 1e-12){
    checkStr = "FAILED";
  }
  std::cout << checkStr << " "
	    << std::setprecision(14)
	    << rNorm1 << " " << trueNorm1 << " "
	    << rNorm2 << " " << trueNorm2 << std::endl;
}

int main(int argc, char *argv[])
{
  // reference state (just for testing)
  fom_state_t yRef(fomSize);
  pressio::ops::fill(yRef, 0.0);

  // basis
  Eigen::MatrixXd phi(fomSize, romSize);
  // phi is just for testing
  phi.setConstant(1.0);
  decoder_t decoderObj(phi);

  testSteadyResidualPolicy(yRef, decoderObj);

  return 0;
}
