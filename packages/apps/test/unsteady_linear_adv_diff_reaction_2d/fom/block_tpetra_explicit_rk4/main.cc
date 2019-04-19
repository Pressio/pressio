
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYLINADVDIFFREACTION2D"
#include "../gold_states_explicit.hpp"

constexpr double eps = 1e-12;
std::string checkStr {"PASSED"};

// template <typename T>
// void checkSol(const T & y,
// 	      const std::vector<double> & trueS){
//   for (size_t i=0; i<trueS.size(); i++){
//     if (std::abs(y[i] - trueS[i]) > eps) checkStr = "FAILED";
//   }
// }

constexpr bool do_print = false;

struct Observer{
  Observer() = default;

  template <typename T>
  void operator()(size_t step,
  		  double t,
  		  const T & y)
  {
    // if (do_print){
    //   if (step % 5 == 0){
    // 	assert( y.localSize() == y.globalSize() );
    // 	std::ofstream file;
    // 	file.open( "sol_" + std::to_string(step) + ".txt" );
    // 	for(auto i=0; i < y.localSize(); i++){
    // 	  file << std::fixed << std::setprecision(14) << y[i] ;
    // 	  file << std::endl;
    // 	}
    // 	file.close();
    //   }
    // }
  }
};

int main(int argc, char *argv[]){
  using app_t		= rompp::apps::UnsteadyLinAdvDiffReac2dBlockTpetra;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

  // scope guard needed (MPI init within trilinos)
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    // create app object
    const int Nx = 7, Ny = 10; //Nx*2-1;
    app_t appobj(Comm, Nx, Ny);
    appobj.setup();
    auto y0n = appobj.getInitialState();
    y0n.putScalar(1);
    auto r0n = appobj.residual(y0n, zero);

    auto r0n_v = r0n.getVectorView();
    auto r0n_v2 = r0n_v.getData();

    for(auto i=0; i<r0n_v2.size(); ++i)
      std::cout << i << " " << std::setprecision(8) << r0n_v2[i] << std::endl;


    // using ode_state_t = rompp::core::Vector<app_state_t>;
    // using ode_res_t   = rompp::core::Vector<app_residual_t>;
    // ode_state_t y(y0n);
    // ode_res_t r(r0n);

    // constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
    // using stepper_t = rompp::ode::ExplicitStepper<
    //   ode_case, ode_state_t, app_t, ode_res_t>;
    // stepper_t stepperObj(y, appobj, r);

    // // integrate in time
    // scalar_t dt = 0.001; scalar_t fint = 0.0010;
    // auto Nsteps = static_cast<unsigned int>(fint/dt);
    // rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps);
    // y.data()->Print(std::cout << std::setprecision(14));
    // {
    //   using namespace rompp::apps::test;
    //   checkSol(y,
    // 	       LinAdvDiffReac2dExpGoldStates<ode_case>::get(Nx, Ny, dt, fint));
    // }
  }

  std::cout << checkStr << std::endl;
  return 0;
}


  // if (do_print){
  //   auto X = appobj.getX(); auto Y = appobj.getY();
  //   std::ofstream file; file.open( "xy.txt" );
  //   for(auto i=0; i < (Nx-2)*Ny; i++){
  //     file << std::fixed << std::setprecision(14)
  // 	   << (*X)[i] << " " << (*Y)[i];
  //     file << std::endl;
  //   }
  //   file.close();
  // }
