
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "../testing_apps/apps.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(int rank, const T & y,
	      const std::vector<double> & trueS)
{
  auto y_v = y.getData();

  int shift = 0;
  if (rank==1) shift = 5;
  else if (rank==2) shift = 10;
  else if (rank==3) shift = 15;

  for (auto i=0; i<5; i++){
    if (std::abs(y_v[i] - trueS[i+shift]) > eps)
      checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using app_t		= pressio::apps::Burgers1dTpetra;
  using scalar_t	= typename app_t::scalar_type;
  using state_t	= typename app_t::state_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    if (Comm->getSize() != 4){
      checkStr = "FAILED";
      return 0;
    }

    std::vector<double> mu({5.0, 0.02, 0.02});
    const int Ncells = 20;
    app_t appObj(mu, Ncells, Comm);
    auto & y0 = appObj.getInitialState();

    state_t y(pressio::ops::clone(y0));
    auto stepperObj = pressio::ode::create_forward_euler_stepper(y,appObj);

    // integrate in time
    scalar_t fint = 35;
    scalar_t dt = 0.01;
    auto Nsteps = static_cast<::pressio::ode::step_count_type>(fint/dt);
    pressio::ode::advance_n_steps(stepperObj, y, 0.0, dt, Nsteps);
    {
      using namespace pressio::apps::test;
      checkSol(rank, y, Burgers1dExpGoldStatesEuler::get(Ncells, dt, fint));
    }
  }

  std::cout << checkStr << std::endl;
  return 0;
}
