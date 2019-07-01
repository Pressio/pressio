
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTION2D"
#include <fstream>

using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReac2dEigen;
using scalar_t		= typename app_t::scalar_type;
using app_state_t	= typename app_t::state_type;
using app_rhs_t	= typename app_t::velocity_type;
using app_jacobian_t	= typename app_t::jacobian_type;

using ode_state_t = rompp::containers::Vector<app_state_t>;
using ode_res_t   = rompp::containers::Vector<app_rhs_t>;
using ode_jac_t   = rompp::containers::Matrix<app_jacobian_t>;

using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
using uint_t		= unsigned int;

constexpr double eps = 1e-8;
constexpr auto t0	= static_cast<scalar_t>(0);

struct FomObserver{
  int colInd_{0};
  eig_dyn_mat A_;

  void resizeRows(int nRows){
    A_.resize(nRows, 0);
  }

  template <typename T>
  void operator()(size_t step,
		  scalar_t t,
		  const T & y){
    A_.conservativeResize(Eigen::NoChange,
			  A_.cols()+1);
    for (auto i=0; i<A_.rows(); i++)
      A_(i,colInd_) = y[i];
    colInd_++;
  }
};


struct FomRunner{
  const int Nx_ = {};
  const int Ny_ = {};
  FomObserver observer_;

  FomRunner(int Nx, int Ny)
    : Nx_{Nx}, Ny_{Ny}, observer_{}{}

  const eig_dyn_mat & getSnapshots() const {
    return observer_.A_;
  }

  ode_state_t run(scalar_t dt, uint_t Nsteps )
  {
    app_t appobj(Nx_, Ny_);
    appobj.setup();
    const auto totDofs = appobj.getUnknownCount();
    const auto y0n = appobj.getInitialState();
    ode_state_t y(y0n);

    constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
    using stepper_t = rompp::ode::ImplicitStepper<
      ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
    stepper_t stepperObj(y, appobj);

    // define solver
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<
      rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
    rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
    solverO.setTolerance(1e-15);
    solverO.setMaxIterations(500);

    // integrate in time
    observer_.resizeRows(totDofs);
    rompp::ode::integrateNSteps(stepperObj, y, t0, dt, Nsteps, observer_, solverO);

    return y;
  }//run
};


void readMatrixFromFile(std::string filename,
			std::vector<std::vector<double>> & A0,
			int ncols){
  assert( A0.empty() );
  std::ifstream source;
  source.open( filename, std::ios_base::in);
  std::string line, colv;
  std::vector<double> tmpv(ncols);
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (int i=0; i<ncols; i++){
      in >> colv;
      tmpv[i] = atof(colv.c_str());
    }
    A0.emplace_back(tmpv);
  }
  source.close();
}


int main(int argc, char *argv[]){

  std::string checkStr {"PASSED"};

  // set parameters to use for runs
  constexpr int Nx = 11, Ny = Nx*2-1;
  constexpr scalar_t dt = 0.1;
  constexpr auto Nsteps = static_cast<unsigned int>(10);
  //constexpr scalar_t fint = Nsteps*dt;

  // run FOM and collect snapshots
  FomRunner fom(Nx, Ny);
  const auto fomY = fom.run(dt, Nsteps);

  // get snapshots
  const auto & S = fom.getSnapshots();
  // do SVD and compute basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  std::cout << std::setprecision(15) << U << std::endl;
  // std::cout << S << std::endl;
  std::cout << "Done with snapshots" << std::endl;

  // std::ofstream file;
  // file.open( "phi.txt" );
  // file << std::fixed
  //      << std::setprecision(15)
  //      << U
  //      << std::endl;
  // file.close();

  // read gold basis from file
  std::vector<std::vector<double>> goldU;
  readMatrixFromFile("gold_basis.txt", goldU, Nsteps+1);

  // check that computed matches gold
  if ( (size_t) goldU.size() != (size_t) U.rows() )
    checkStr = "fail";

  std::cout.precision(15);
  for (auto i=0; i<U.rows(); i++)
    for (auto j=0; j<U.cols(); j++){
      auto err = std::abs(goldU[i][j] - U(i,j));
      std::cout << "gold = " << goldU[i][j]
		<< " computed = " << U(i,j)
		<< " |err| = " << err
		<< std::endl;

      if ( err > eps )
	checkStr = "FAILED";
    }

  std::cout << checkStr << std::endl;
  return 0;
}
