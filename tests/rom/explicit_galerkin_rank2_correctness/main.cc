
#include "pressio_rom_galerkin.hpp"

struct MyFakeApp
{
  using scalar_type	  = double;
  using state_type	  = Eigen::Matrix<scalar_type, -1,-1, Eigen::ColMajor>;
  using velocity_type     = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

  int nFields_ = {};
  int N_ = {};

  MyFakeApp(int N, int nFields)
    : nFields_(nFields), N_(N){}

  velocity_type createVelocity() const
  {
    velocity_type f(N_, nFields_);
    return f;
  }

  void velocity(const state_type & u,
  		const scalar_type t,
		velocity_type & f) const
  {
    std::cout << t << std::endl;
    for (auto j=0; j<f.cols(); ++j){
      for (auto i=0; i<f.rows(); ++i){
	f(i,j) = j==0 ? 1.+t : j==1 ? 2.+t : 3.+t;
      }
    }
  }
};

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});
  std::string checkStr {"PASSED"};

  using fom_t	= MyFakeApp;
  using scalar_t = typename fom_t::scalar_type;
  using native_state_t	= typename fom_t::state_type;
  using native_dmat_t	= typename fom_t::dense_matrix_type;

  using fom_state_t = pressio::containers::MultiVector<native_state_t>;
  using basis_t	= pressio::containers::experimental::MultiVectorSet<native_dmat_t>;
  using decoder_t = pressio::rom::LinearDecoder<basis_t, fom_state_t>;
  using rom_state_t = pressio::containers::MultiVector<Eigen::MatrixXd>;

  const int gridSize  = 8;
  // e.g. density, velo, temperature
  const int numFields = 3;
  // number of modes for each physical quantity
  const int numModes  = 4;

  // create fom object
  fom_t appobj(gridSize,numFields);

  // create the basis
  basis_t phi(3, gridSize, numModes);
  scalar_t vals[3] = {1., 2., 3.};
  for (int i=0; i<numFields; ++i){
    ::pressio::ops::fill(phi(i), vals[i] );
  }

  // decoder object
  decoder_t decoder(phi);

  // reference state
  native_state_t yRef(gridSize, numFields);
  yRef.setConstant(0.);

  // rom state
  rom_state_t romState(numModes, numFields);
  ::pressio::ops::fill(romState, 0.0);

  // create problem
#ifdef EULER
  using ode_tag = pressio::ode::explicitmethods::Euler;
#endif
#ifdef RK4
  using ode_tag = pressio::ode::explicitmethods::RungeKutta4;
#endif
  using pressio::rom::galerkin::createDefaultProblem;
  auto Problem = createDefaultProblem<ode_tag>(appobj, decoder, romState, yRef);

  pressio::rom::galerkin::solveNSteps(Problem, romState, 0.0, 1., 1);

  /*
    phi[:,:,0] = all 1s
    phi[:,:,1] = 2s
    phi[:,:,2] = 3s

    //------------------
    for forward Euler:
    //------------------
       f[:,0] f[:,1] f[:,2]
 	1    2    3
	1    2    3
	1    2    3
	1    2    3
	1    2    3
	1    2    3
	1    2    3
	1    2    3

       B= phi^T f =  (where we consider one slice of phi per field)
		[8 32 72]
		[8 32 72]
		[8 32 72]
		[8 32 72]

       xNew[:,0] = xOld[:,0] + dt * B[:,0]
		 = [0 0 0 0] + 1. * [8 8 8 8]
       xNew[:,1] = xOld[:,1] + dt * B[:,1]
                 = [0 0 0 0] + 1. * [32 32 32 32]
       xNew[:,1] = xOld[:,1] + dt * B[:,2]
	         = [0 0 0 0] + 1. * [72 72 72 72]

    //------------------
    for rk4:
    //------------------
    xNew[:,0] = xOld[:,0]
	+ (1/6)*dt*(
  	phi[:,:,0]^T*k1[:,0] +
	phi[:,:,0]^T*2*k2[:,0] +
	phi[:,:,0]^T*2*k3[:,0] +
	phi[:,:,0]^T*k4[:,0]
	)

    k1
      1    2    3
      1    2    3
      ...

    k2
     1+dt/2  2+dt/2  3+dt/2
     1+dt/2  2+dt/2  3+dt/2
      ...

    k3
     1+dt/2  2+dt/2  3+dt/2
     1+dt/2  2+dt/2  3+dt/2
      ...

    k4
      2    3    4
      2    3    4
      ...

    xNew[:,0] = [0 0 0 0 ]
	+ 1/6 * 1 * ([8 8 8 8]
		   + [24 24 24 24]
		   + [24 24 24 24]
		   + [16 16 16 16])
              = [12 12 12 12]

    xNew[:,1] = [40 40 40 40]
    xNew[:,2] = [84 84 84 84]
   */

#ifdef EULER
  Eigen::MatrixXd goldRom(numModes, numFields);
  goldRom << 8.,32.,72.,8.,32.,72.,8.,32.,72.,8.,32.,72.;
#endif
#ifdef RK4
  Eigen::MatrixXd goldRom(numModes, numFields);
  goldRom << 12.,40.,84.,12.,40.,84.,12.,40.,84.,12.,40.,84.;
#endif

  std::cout << *romState.data() << std::endl;
  std::cout << goldRom << std::endl;
  if (!romState.data()->isApprox(goldRom)) checkStr = "FAILED";

  std::cout << checkStr <<  std::endl;
  return 0;
}
