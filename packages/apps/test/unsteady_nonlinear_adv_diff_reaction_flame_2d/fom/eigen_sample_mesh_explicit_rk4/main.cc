
#include "CORE_ALL"
#include "ODE_ALL"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../gold_states_explicit.hpp"

constexpr double eps = 1e-10;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y,
	      const std::vector<double> & trueS){
  if (trueS.empty()) {
    std::cout << " true solution not found, empty " << std::endl;
    checkStr = "FAILED";
  }
  for (size_t i=0; i<trueS.size(); i++){
    const auto err = std::abs(y[i] - trueS[i]);
    std::cout << std::fixed << std::setprecision(15)
	      << " true = " << trueS[i]
	      << " y = " << y[i]
	      << " err = " << err
	      << std::endl;
    if ( err > eps or std::isnan(y[i])) checkStr = "FAILED";
  }
}

constexpr bool do_print = true;

struct Observer{
  Observer() = default;

  template <typename T>
  void operator()(size_t step,
  		  double t,
  		  const T & y)
  {
    if (do_print){
      if (step % 50 == 0){
    	std::ofstream file;
    	file.open( "sol_" + std::to_string(step) + ".txt" );
    	for(auto i=0; i < y.size(); i++){
    	  file << std::fixed << std::setprecision(14) << y[i] ;
    	  file << std::endl;
    	}
    	file.close();
      }
    }
  }
};


template <typename T>
void readMeshGraph(std::string filename,
		   std::vector<T> & A){
  std::ifstream file(filename);
  std::string line;

  T lineGIDs;
  while (getline(file, line)){
    std::istringstream ss(line);

    int entry;
    int j = 0;
    while (ss >> entry){
      lineGIDs[j] = entry;
      j++;
    }
    A.emplace_back(lineGIDs);
  }

  for (auto & it : A){
    for (auto & it2 : it){
      std::cout << " " << it2;
    }
    std::cout << std::endl;
  }
}

void readMappingGIDs(std::string filename,
		     std::vector<std::array<int,2>> & A){
  std::ifstream file(filename);
  std::string line;
  std::array<int,2> lineGIDs;
  while (getline(file, line)){
    std::istringstream ss(line);

    int entry;
    int j = 0;
    while (ss >> entry){
      lineGIDs[j] = entry;
      j++;
    }
    A.emplace_back(lineGIDs);
  }

  for (auto & it : A){
    for (auto & it2 : it){
      std::cout << " " << it2;
    }
    std::cout << std::endl;
  }
}



struct updateOps{
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
  using scalar_t	= typename app_t::scalar_type;
  using s_t		= typename app_t::state_type;
  using r_t		= typename app_t::residual_type;

  static void do_update(s_t & v, const scalar_t c,
  			const r_t & v0, const scalar_t a,
  			const r_t & v1, const scalar_t b){
    for (auto i=0; i<v.size(); ++i)
      v[i] = c*v[i] + a*v0[i] + b*v1[i];
  }

  static void do_update(s_t & v,
  			const r_t & v0, const scalar_t a,
  			const r_t & v1, const scalar_t b){
    for (auto i=0; i<v.size(); ++i)
      v[i] = a*v0[i] + b*v1[i];
  }

  static void do_update(s_t & v,
			const r_t & v1, const scalar_t & b,
			const r_t & v2, const scalar_t & c,
			const r_t & v3, const scalar_t & d,
			const r_t & v4, const scalar_t & e) {
    for (auto i=0; i<v.size(); ++i)
      v[i] = b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }

  static void do_update(s_t & v, const scalar_t & a,
			const r_t & v1, const scalar_t & b,
			const r_t & v2, const scalar_t & c,
			const r_t & v3, const scalar_t & d,
			const r_t & v4, const scalar_t & e) {
    for (auto i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }
};

struct myops{
  using update_op = updateOps;
};


int main(int argc, char *argv[]){
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
  using scalar_t	= typename app_t::scalar_type;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

  /*
    mesh.dat: contains list of gid and for each gid the neighbors needed
	to compute the residual at that location.
	* the neighbors are ordered as east, north, west, south.
	* if a neighbor is = -1, it means the current grid point is near boundary
	and since neighbors are ordered, we know which boundary is near this point.
	* if a grid point is meant to only hold state, then all neighbors are = -2

    gappy_to_full_gid_map.dat: it is a 2-column file where the first column
	contains the gids of each grid point enumerated according to the sample mesh
	while the second column contains the corresponding gid of each grid within
	the full mesh.
   */

  typename app_t::graph_t meshGraph;
  readMeshGraph("mesh.dat", meshGraph);

  typename app_t::gids_map_t gToFgidMap;
  readMappingGIDs("gappy_to_full_gid_map.dat", gToFgidMap);

  constexpr int Nx = 40, Ny = 20;
  app_t appobj(Nx, Ny, meshGraph, gToFgidMap);
  appobj.setup();
  const auto y0n = appobj.getInitialState();
  const auto r0n = appobj.residual(y0n, zero);

  std::cout << r0n << std::endl;

  // if (do_print){
    auto X = appobj.getX(); auto Y = appobj.getY();
    std::ofstream file; file.open( "xy.txt" );
    for(auto i=0; i < X.size(); i++){
      std::cout << std::fixed
		<< std::setprecision(15)
		<< X[i] << " " << Y[i] << " "
		<< y0n[i] << " "
		<< r0n[i*4] << " " << r0n[i*4+1] << " "
		<< r0n[i*4+2] << " " << r0n[i*4+3];
      std::cout << std::endl;

      file << std::setprecision(14)
  	   << X[i] << " " << Y[i];
      file << std::endl;
    }
    file.close();
  // }

  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_residual_t>;
  ode_state_t y(y0n);
  ode_res_t r(r0n);

  constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
  using stepper_t = rompp::ode::ExplicitStepper<
    ode_case, ode_state_t, app_t, ode_res_t, scalar_t, myops>;
  stepper_t stepperObj(y, appobj, r);

  // integrate in time
  constexpr scalar_t dt = 0.00001;
  constexpr auto Nsteps = 1000;
  constexpr scalar_t fint = Nsteps*dt;
  Observer obs;
  rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs);
  //std::cout << std::fixed << std::setprecision(15) << *y.data() << std::endl;

  for(auto i=0; i < X.size(); i++){
    std::cout << std::fixed
  	      << std::setprecision(15)
  	      << X[i] << " " << Y[i] << " "
  	      << (*y.data())[i*4] << " "
  	      << (*y.data())[i*4+1] << " "
  	      << (*y.data())[i*4+2] << " "
  	      << (*y.data())[i*4+3];
    std::cout << std::endl;
  }


  // {
  //   using namespace rompp::apps::test;
  //   checkSol(y,
  // 	     NonLinAdvDiffReacFlame2dExpGoldStates<ode_case>::get(Nx,
  // 								  Ny,
  // 								  dt,
  // 								  fint));
  // }

  std::cout << checkStr << std::endl;
  return 0;
}
