
#include "pressio_rom_lspg.hpp"

struct MyMasker
{
  using scalar_type = double;

  int N_ = {};
  MyMasker(int N) : N_(N){}

  Eigen::VectorXd createApplyMaskResult(const Eigen::VectorXd & operand) const
  {
    return Eigen::VectorXd(N_);
  }
  Eigen::MatrixXd createApplyMaskResult(const Eigen::MatrixXd & operand) const
  {
    Eigen::MatrixXd A;
    return A;
  }

  void applyMask(const Eigen::VectorXd & operand,
     const scalar_type & t,
     Eigen::VectorXd & result) const
  {
    for (int i=0; i<N_; ++i) result(i) = operand(i) + 1.5;
  }

  void applyMask(const Eigen::MatrixXd &,
     const scalar_type & t,
     Eigen::MatrixXd &) const{}
};

struct MyPreconditioner
{
  using scalar_type = double;

  void applyPreconditioner(const Eigen::VectorXd & yState,
         const scalar_type & time,
         Eigen::VectorXd & operand) const
  {
    for (int i=0; i<operand.size(); ++i) operand(i) += 0.5;
  }

  void applyPreconditioner(const Eigen::VectorXd & yState,
         const scalar_type & time,
         Eigen::MatrixXd & C) const
  {}
};

struct MyFakeApp
{
  int N_ = {};
public:
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  MyFakeApp(int N) : N_(N){}

  velocity_type createVelocity() const{ return state_type(N_);}
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    dense_matrix_type A;
    return A;
  }

  void velocity(const state_type & state,
		const double & time,
		velocity_type & f) const
  {
    for (int i=0; i<N_; ++i) f(i) = 1;
  }
  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const
  {}
};

struct MyFakeSolver
{
  MyFakeSolver(){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state){}
};


constexpr int fomSize = 10;
constexpr int romSize = 4;

using sc_t		= double;
using native_state_t	= Eigen::VectorXd;
using native_dmat_t	= Eigen::MatrixXd;
using fom_state_t	= pressio::containers::Vector<native_state_t>;
using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
using decoder_t		= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

using lspg_state_t    = pressio::containers::Vector<Eigen::VectorXd>;
using lspg_residual_t = fom_state_t;

using preconditioner_t = MyPreconditioner;
using masker_t = MyMasker;


void testUnsteadyResidualPolicy(fom_state_t & yRef,
				decoder_t & decoderObj)
{
  // create fom object
  MyFakeApp fomObj(fomSize);

  std::string checkStr {"PASSED"};

  // we test here for BDF1
  constexpr std::size_t numAuxStates = 1;
  using ode_tag  = pressio::ode::implicitmethods::BDF1;

  // ------------------------------------
  // construct what we need to construct
  // the residual policy object
  // ------------------------------------
  using fom_recon = pressio::rom::FomStateReconstructor<sc_t,fom_state_t,decoder_t>;
  fom_recon fomReconstructor(yRef, decoderObj);

  using fom_states_manager_t = pressio::rom::ManagerFomStates<
    ::pressio::rom::UnsteadyImplicit, fom_state_t, fom_recon, void, numAuxStates+1>;
  fom_states_manager_t fomStatesMngr(fomReconstructor, yRef);

  using lspg_residual_policy =
    ::pressio::rom::lspg::impl::unsteady::ResidualPolicyContinuousTimeApi<
      fom_state_t, fom_states_manager_t, void>;

  using lspg_precond_residual_policy =
    ::pressio::rom::lspg::decorator::Preconditioned<preconditioner_t, lspg_residual_policy>;

  using lspg_masked_residual_policy =
    ::pressio::rom::lspg::decorator::Masked<masker_t, lspg_residual_policy>;

  const lspg_residual_policy rPol(fomStatesMngr);

  preconditioner_t Preconditioner;
  const lspg_precond_residual_policy rPolPrecond(Preconditioner, fomStatesMngr);

  masker_t Masker(fomSize);
  const lspg_masked_residual_policy rPolMasked(Masker, fomObj, fomStatesMngr);

  // ------------------------------------
  // ***** here we do the test *****
  // ------------------------------------
  using lspg_aux_states = pressio::ode::implicitmethods::StencilStatesManager<lspg_state_t, numAuxStates>;

  lspg_state_t romState(romSize);
  lspg_aux_states romAuxStates(romState);

  lspg_residual_t residual(fomSize);

  // set romY_n+1 = {1,2,3,4}
  romState(0) = 1.;
  romState(1) = 2.;
  romState(2) = 3.;
  romState(3) = 4.;

  // set romY_n = {4,1,2,1}
  auto &  romYn = romAuxStates(pressio::ode::n());
  romYn(0) = 4.;
  romYn(1) = 1.;
  romYn(2) = 2.;
  romYn(3) = 1.;

  /* the compute call below should first compute the residual for BDF1:
	R = phi*romY_n - phi*romY_n-1 - dt*f( phi*romY_n)
     and then call:
	applyPreconditioned(R) on the fomObj

     where we set f to be all 1, dt = .5, phi is all 1.
     - the apply precon takes R and mimics a modification by adding 0.5 to each entry
     - the apply mask takes R and mimics a modification by adding 1.5 to each entry
     here we don't care if this is meaningful, we just want to test things.

     so we should have:
     fomY_n = [1 1 1 1|   1           10
	      |1 1 1 1|   2           10
	      |1 1 1 1|   3     =     10
	      |1 1 1 1|   4           10
	      |...	              ...
	      |1 1 1 1]		      10

     fomY_n-1 = [1 1 1 1|   4           8
		|1 1 1 1|   1           8
		|1 1 1 1|   2     =     8
		|1 1 1 1|   1           8
		|...	                ...
		|1 1 1 1]  	        8

     so before preconditioning, we have
	R  = [1.5 1.5 ... 1.5]

     the precond should give:
	R  = [2 2 ... 2]

     the mask should give:
	R  = [3 3 ... 3]
  */

  sc_t rNorm1 = {};
  sc_t rNorm2 = {};
  sc_t rNorm3 = {};

  rPol.compute<ode_tag>(romState, romAuxStates, fomObj,
    0.0, 0.5, 1,residual);
  rNorm1 = residual.data()->norm();

  rPolPrecond.compute<ode_tag>(romState,romAuxStates,fomObj,
    0.0, 0.5, 1,residual);
  rNorm2 = residual.data()->norm();

  rPolMasked.compute<ode_tag>(romState, romAuxStates, fomObj,
    0.0, 0.5, 1,residual);
  rNorm3 = residual.data()->norm();

  const auto trueNorm1 = std::sqrt( (1.5*1.5)*(double)fomSize );
  const auto trueNorm2 = std::sqrt( (2.0*2.0)*(double)fomSize );
  const auto trueNorm3 = std::sqrt( (3.0*3.0)*(double)fomSize );

  const auto err1 = std::abs( rNorm1 - trueNorm1 );
  const auto err2 = std::abs( rNorm2 - trueNorm2 );
  const auto err3 = std::abs( rNorm3 - trueNorm3 );

  if ( err1 > 1e-12 or err2 > 1e-12 or err3>1e-12){
    checkStr = "FAILED";
  }
  std::cout << checkStr << " "
	    << std::setprecision(14)
	    << rNorm1 << " " << trueNorm1 << " "
	    << rNorm2 << " " << trueNorm2 << " "
	    << rNorm3 << " " << trueNorm3 <<  std::endl;
}

int main(int argc, char *argv[])
{
  fom_state_t yRef(fomSize);
  pressio::ops::fill(yRef, 0.0);

  // basis
  Eigen::MatrixXd phi(fomSize, romSize);
  // phi is just for testing
  phi.setConstant(1.0);
  decoder_t decoderObj(phi);

  testUnsteadyResidualPolicy(yRef, decoderObj);

  return 0;
}
