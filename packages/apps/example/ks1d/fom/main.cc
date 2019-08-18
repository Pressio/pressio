
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "APPS_KS1D"

#include <iostream>
#include <fstream>

constexpr double eps = 1e-12;

template <typename state_t>
struct observer{
  using matrix_t = Eigen::MatrixXd;

  size_t state_size_ {};
  matrix_t A_;
  matrix_t J_;
  size_t count_ {};
  state_t y0_;
  state_t yIncr_;
  state_t Jbar_;
  int Nsamp_;

  observer(int N, int state_size, const state_t & y0)
    : state_size_(state_size),
      A_(state_size, N+1), //+1 to store also init cond
	  J_(4, N+1),
      y0_(y0),
      yIncr_(state_size),
	  Jbar_(4),
	  Nsamp_(N){}

  void operator()(size_t step,
  		  double t,
  		  const state_t & y){
    yIncr_ = y; // - y0_; // Just want to plot state for now
    this->storeInColumn(yIncr_, count_);
    this->computeQoI(y);

    count_++;
  }

  void storeInColumn(const state_t & y, int j){
    for (auto i=0; i<y.size(); i++)
      A_(i,j) = y(i);
  }

  void computeQoI(const state_t & y)
  {
	double qoi_array[4];

	// point quantities
	qoi_array[0] = y[state_size_/4];
	qoi_array[1] = y[state_size_/4] * y[state_size_/4];

	// Spatial averages
	qoi_array[2] = 0.0;
	qoi_array[3] = 0.0;
	for (size_t i=0; i < state_size_; i++)
	{
	  qoi_array[2] += y[i] / (state_size_ + 2); // +2 for boundary nodes
	  qoi_array[3] += y[i] * y[i] / (state_size_ + 2);
	}

	for (int i=0; i < 4; i++)
	{
	  // save quantities of interest
	  J_(i, count_) = qoi_array[i];
	  // accumulate time average quantity of interest
	  Jbar_(i) += (1.0/(Nsamp_+1.0)) * qoi_array[i];
	}
  }

  void printAll() const
  {
	  std::ofstream myfile;
	  myfile.open("fom.dat");
	  myfile << std::setprecision(14) << A_ << std::endl;
	  myfile.close();

	  std::ofstream myQoIfile;
	  myQoIfile.open("QoI.dat");
	  myQoIfile << std::setprecision(14) << J_ << std::endl;
	  myQoIfile.close();

  }

  void printFinal() const
  {
	  std::ofstream myfile;
	  myfile.open("fom_final.dat");

	  for (auto i=0; i<y0_.size(); i++)
	  {
		  myfile << std::setprecision(14) << A_(i,count_-1) << std::endl;
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
void read_inputs(scalar_t & c, scalar_t & dt , scalar_t & T, scalar_t & L, int & nx )
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
				// Do nothing
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
  using app_t		= pressio::apps::KS1dEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  //-------------------------------
  // Read inputs
  scalar_t c = 0.0;
  scalar_t dt = -1.0;
  scalar_t fint = -1.0;
  scalar_t L = -1.0;
  int nx;

  read_inputs( c,  dt , fint, L, nx );

  // create app object
  int Nnodes = nx;
  Eigen::Vector3d mu(c, L, 0.0);
  app_t appObj(mu, Nnodes);
  appObj.setup();
  auto & y0n = appObj.getInitialState();

  // types for ode
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::Matrix<app_jacob_t>;

  ode_state_t y(y0n);

  // define auxiliary stepper
  using aux_stepper_t = pressio::ode::ImplicitStepper<
    pressio::ode::ImplicitEnum::Euler,
    ode_state_t, ode_res_t, ode_jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  // nonimal stepper
  constexpr auto ode_case = pressio::ode::ImplicitEnum::BDF2;
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

  // define solver
  using lin_solver_t = pressio::solvers::iterative::EigenIterative<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  pressio::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
  solverO.setTolerance(1e-13);
  solverO.setMaxIterations(200);

  // integrate in time
//  scalar_t fint = 1000.0;
//  scalar_t dt = 0.1;
  auto Nsteps = static_cast<unsigned int>(fint/dt);

  // define observer
  observer<ode_state_t> Obs(Nsteps, Nnodes, y);

  pressio::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, Obs, solverO);
  Obs.printAll(); // print all states and QoIs
  Obs.printFinal(); // print final state and time-averaged QoIs

  return 0;
}
