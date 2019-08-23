
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_KS1D"
#include "utils_eigen.hpp"

// constexpr double eps = 1e-12;

template <typename state_t,typename lspg_t>
struct observer{
  using matrix_t = Eigen::MatrixXd;

  size_t state_size_ {};
  matrix_t A_;
  matrix_t J_;
  size_t count_ {};
  state_t y0_;
  lspg_t lspg_;
  state_t Jbar_;
  int Nsamp_;

  observer(int N, int state_size, const state_t & y0, const lspg_t & lspg)
    : state_size_(state_size),
      A_(state_size, N+1), //+1 to store also init cond
	  J_(4, N+1),
      y0_(y0),
	  lspg_(lspg),
      Jbar_(4),
      Nsamp_(N){}

  void operator()(size_t step,
  		  double t,
  		  const state_t & y){
    this->storeInColumn(y, count_);
    count_++;
  }

  void storeInColumn(const state_t & y, int j){
    for (auto i=0; i<y.size(); i++)
      A_(i,j) = y(i);
  }

  //TODO project to full order, then compute QoIs
  void computeQoIs(){

	std::cout << "Computing QoIs" << std::endl;

	state_t yROM(state_size_);

	for (size_t j=0; j < Nsamp_+1; j++) {
		// Get gen coord
        for (size_t i=0; i<state_size_; i++ ) {
        	  yROM(i) = A_(i,j);
        }

		// Reconstruct FOM state
		auto y = lspg_.yFomReconstructor_(yROM);

		size_t fom_size = y.size();

		double qoi_array[4];

		// point quantities
		qoi_array[0] = y[fom_size/4];
		qoi_array[1] = y[fom_size/4] * y[fom_size/4];

		// Spatial averages
		qoi_array[2] = 0.0;
		qoi_array[3] = 0.0;
		for (size_t i=0; i < fom_size; i++)
		{
		  qoi_array[2] += y[i] / (fom_size + 2); // +2 for boundary nodes
		  qoi_array[3] += y[i] * y[i] / (fom_size + 2);
		}

		for (int i=0; i < 4; i++)
		{
		  // save quantities of interest
		  J_(i, j) = qoi_array[i];
		  // accumulate time average quantity of interest
		  Jbar_(i) += (1.0/(Nsamp_+1.0)) * qoi_array[i];
		}

	}


  }

  void printAll() const
  {
	  std::ofstream myfile;
	  myfile.open("gen_coords.dat");
	  myfile << std::setprecision(14)  << A_ << std::endl;
	  myfile.close();

	  std::ofstream myQoIfile;
	  myQoIfile.open("QoI.dat");
	  myQoIfile << std::setprecision(14) << J_ << std::endl;
	  myQoIfile.close();


  }

  void printFinal() const
  {

	  // compute the fom corresponding to our rom final state
	  state_t yROM(y0_);
	  for (auto i=0; i<y0_.size(); i++)
	  {
		yROM(i) = A_(i,count_-1);
	  }

	  auto yFomFinal = lspg_.yFomReconstructor_(yROM);

	  std::ofstream myfile;
	  myfile.open("rom_final.dat");

	  for (auto i=0; i<yFomFinal.size(); i++)
	  {
		  myfile << std::setprecision(14) << yFomFinal(i) << std::endl;
	  }

	  myfile.close();

	  std::ofstream myQoIfile;
	  myQoIfile.open("QoI_avg.dat");
	  for (auto i=0; i<Jbar_.size(); i++)
	  {
	      myQoIfile << std::setprecision(14) << Jbar_(i) << std::endl;
	  }
	  myQoIfile.close();
  }


};

template <typename scalar_t>
void read_inputs(scalar_t & c, scalar_t & nu, scalar_t & dt , scalar_t & T, scalar_t & L, int & nx , int & romSize )
{
	// Read in parameters and grid size


	std::cout << "Reading Inputs" << std::endl;

	std::ifstream infile("inputs.txt");

	std::string line;


	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string name;
		scalar_t num;
		if (!(iss >> name >> num)) { break; } // error

		std::cout << name << " " << num << std::endl;

		try
		{
			if (name.compare("c:")==0)
			{
				c = num;
			}
			else if (name.compare("nu:")==0)
			{
				nu = num;
			}
			else if (name.compare("T:")==0)
			{
				T = num;
			}
			else if (name.compare("dt:")==0)
			{
				dt = num;
			}
			else if (name.compare("L:")==0)
			{
				L = num;
			}
			else if (name.compare("nx:")==0)
			{
				nx = (int) num;
			}
			else if (name.compare("romSize:")==0)
			{
				romSize = (int) num;
			}
			else
			{
				throw 20;
			}
		}
		catch (int e)
		{
			std::cout << "Do not recognize the input name " << name << std::endl;
		}

	}

	infile.close();

}


int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::KS1dEigen;
  using scalar_t	= typename fom_t::scalar_type;

  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using lspg_state_t	= pressio::containers::Vector<eig_dyn_vec>;

  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  //-------------------------------

  // Read inputs
  scalar_t c = 0.0;
  scalar_t nu = 1.0;
  scalar_t dt = -1.0;
  scalar_t fint = -1.0;
  scalar_t L = -1.0;
  int nx;
  int romSize;

  read_inputs( c, nu, dt, fint, L, nx, romSize );

  // create app object
  int numNode = nx;
  Eigen::Vector3d mu(c, nu, L);
  fom_t appobj( mu, numNode);
  appobj.setup();
  auto t0 = static_cast<scalar_t>(0);

  // read from file the jacobian of the decoder

  // store modes computed before from file
  decoder_jac_t phi =
    pressio::apps::test::eigen::readBasis("basis.txt", romSize, numNode);
  const int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // for this problem, my reference state = initial state
  auto & yRef = appobj.getInitialState();

  // define ROM state
  lspg_state_t yROM(romSize);
  // initialize to zero (this has to be done)
  yROM.putScalar(0.0);

  // define LSPG type
  constexpr auto ode_case  = pressio::ode::ImplicitEnum::BDF2;
  using lspg_problem_types = pressio::rom::DefaultLSPGTypeGenerator<
    fom_t, ode_case, decoder_t, lspg_state_t>;
  using lspg_problem_gen_t = pressio::rom::LSPGUnsteadyProblemGenerator<
		  lspg_problem_types>;
  lspg_problem_gen_t lspgProblem(
      appobj, yRef, decoderObj, yROM, t0);

  using lspg_stepper_t = typename lspg_problem_types::lspg_stepper_t;

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // GaussNewton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using gnsolver_t   = pressio::solvers::iterative::GaussNewton<
    lspg_stepper_t, linear_solver_t>;
  gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
  solver.setTolerance(1e-13);
  solver.setMaxIterations(200);

  // integrate in time
  auto Nsteps = static_cast<unsigned int>(fint/dt);

  // define observer
  observer<lspg_state_t,lspg_problem_gen_t> Obs(Nsteps, romSize, yROM, lspgProblem);

  pressio::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, 0.0, dt, Nsteps, Obs, solver);

  Obs.computeQoIs();
  Obs.printAll();
  Obs.printFinal();


  //TODO compute and print QoI(s)


  return 0;
}
