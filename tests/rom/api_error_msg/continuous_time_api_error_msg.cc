
#include "pressio_rom.hpp"

struct ValidApp{
  const int32_t numDof_ = 15;

  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using velocity_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;
  using dense_matrix_type = Eigen::MatrixXd;

public:
    // if FAILMODE1 is defined, exclude this method so to trigger failure
#if defined FAILMODE1
#else
  void velocity(const state_type &,
		const scalar_type & time,
		velocity_type &) const
  {}
#endif

#if defined FAILMODE2
#else
  void applyJacobian(const state_type &,
		     const dense_matrix_type &,
		     const scalar_type &,
		     dense_matrix_type &) const
  {}
#endif

#if defined FAILMODE3
#else
  velocity_type createVelocity() const
  {
    return velocity_type(numDof_);
  }
#endif

#if defined FAILMODE4
#else
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type & B) const{
    return dense_matrix_type(numDof_, 3);
  }
#endif
};

int main(int argc, char *argv[])
{
  using namespace pressio;
  using app_t    = ValidApp;
  using scalar_t	= typename app_t::scalar_type;
  using native_state_t  = typename app_t::state_type;
  using fom_state_t  = pressio::containers::Vector<native_state_t>;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using decoder_jac_t
     = pressio::containers::MultiVector<typename app_t::dense_matrix_type>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;

  app_t appobj;
  decoder_jac_t phi(Eigen::MatrixXd(appobj.numDof_,3));
  decoder_t decoderObj(phi);
  typename app_t::state_type yRef;
  lspg_state_t yROM;

  using ode_name_t = pressio::ode::implicitmethods::Euler;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_name_t, app_t, lspg_state_t, decoder_t>::type;

  lspg_problem lspgProblem(appobj, yRef, decoderObj, yROM);

  return 0;
}
