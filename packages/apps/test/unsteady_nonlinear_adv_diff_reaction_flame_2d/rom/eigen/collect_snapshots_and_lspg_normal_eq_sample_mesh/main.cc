
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../../../fom/gold_states_implicit.hpp"

using scalar_t		= double;
using uint_t		= unsigned int;
using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
constexpr auto zero	= ::rompp::core::constants::zero<scalar_t>();
constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
constexpr auto t0	= zero;

struct FomObserver{
  int snapId_{0};
  eig_dyn_mat A_;
  eig_dyn_vec y0_;

  void resizeRows(int nRows){
    // num of rows = num of dofs
    A_.resize(nRows, 0);
    y0_.resize(nRows,1);
  }

  template <typename T>
  void operator()(size_t step,
		  scalar_t t,
		  const T & y){
    if (step == 0){
      for (auto i=0; i<y0_.rows(); i++)
	y0_[i] = y[i];
    }
    else{
      // snapshots are stored after subtracting init cond
      A_.conservativeResize(Eigen::NoChange,
			    A_.cols()+1);
      for (auto i=0; i<A_.rows(); i++)
	A_(i,snapId_) = y[i]-y0_[i];
      snapId_++;
    }
  }
};


struct FomRunner{
  using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dEigen;
  using app_state_t	= typename app_t::state_type;
  using app_residual_t	= typename app_t::residual_type;
  using app_jacobian_t	= typename app_t::jacobian_type;
  using ode_state_t = rompp::core::Vector<app_state_t>;
  using ode_res_t   = rompp::core::Vector<app_residual_t>;
  using ode_jac_t   = rompp::core::Matrix<app_jacobian_t>;

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

    using stepper_t = rompp::ode::ImplicitStepper<
      ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
    stepper_t stepperObj(y, appobj);

    // define solver
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<
      rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
    rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
    solverO.setTolerance(1e-6);
    solverO.setMaxIterations(100);

    // allocate space in observer to store snapshots
    observer_.resizeRows(totDofs);

    // integrate in time
    rompp::ode::integrateNSteps(stepperObj, y, t0, dt,
				Nsteps, observer_, solverO);

    return y;
  }//run
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

  // for (auto & it : A){
  //   for (auto & it2 : it){
  //     std::cout << " " << it2;
  //   }
  //   std::cout << std::endl;
  // }
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

  // for (auto & it : A){
  //   for (auto & it2 : it){
  //     std::cout << " " << it2;
  //   }
  //   std::cout << std::endl;
  // }
}



struct time_discrete_ops{
  using app_sm_t   = rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
  using graph_t	   = typename app_sm_t::graph_t;
  using gids_map_t = typename app_sm_t::gids_map_t;
  using state_t	   = rompp::Vector<typename app_sm_t::state_type>;
  using resid_t	   = rompp::Vector<typename app_sm_t::state_type>;

  const graph_t & graph_;
  const gids_map_t & gidsMap_;

  time_discrete_ops(const graph_t & graph,
		    const gids_map_t & gidsMap)
    : graph_{graph}, gidsMap_{gidsMap}{
  }

  void time_discrete_residual(const state_t & yn,
			      const std::array<state_t,n> & ynm,
			      fom_resid_t & R,
			      scalar_type dt){
    // // On input: R contains the application RHS, i.e. if
    // // dudt = f(x,u,...), R contains f(...)

    // here we need to use the graph to do proper ops since state and residual are different
  }
}




int main(int argc, char *argv[]){
  std::string checkStr {"PASSED"};

  /*
    run the FOM with full mesh and create basis
  */
  // set parameters to use for runs
  constexpr int Nx = 20, Ny = 10;
  constexpr scalar_t dt = 0.001;
  constexpr auto Nsteps = static_cast<uint_t>(10);

  // run FOM and get snapshots using full mesh case
  // run FOM and collect snapshots
  FomRunner fom(Nx, Ny);
  const auto fomY = fom.run(dt, Nsteps);

  // get snapshots
  const auto & S = fom.getSnapshots();
  std::cout << S << std::endl;

  // compute SVD to create basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  std::cout << std::setprecision(15) << U << std::endl;
  // std::cout << S << std::endl;
  std::cout << "Done with snapshots" << std::endl;

  /*
    run ROM with sample mesh
  */
  {
    // this is the size of the rom
    constexpr int romSize = Nsteps;

    // the type of the app with sample mesh
    using app_sm_t	 = rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;

    // read sample mesh
    typename app_sm_t::graph_t meshGraph;
    readMeshGraph("mesh.dat", meshGraph);
    typename app_sm_t::gids_map_t smGidsToFGidsMap;
    readMappingGIDs("gappy_to_full_gid_map.dat", smGidsToFGidsMap);

    // create app object
    app_sm_t appobj(Nx, Ny, meshGraph, smGidsToFGidsMap);
    appobj.setup();

    // typedefs used for rom
    using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::core::MultiVector<eig_dyn_mat>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // use sample mesh to extract only taget rows from basis
    // ...

    // create decoder
    decoder_jac_t phi(U);
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();

    // define ROM state and set to zero
    lspg_state_t yROM(romSize);
    yROM.putScalar(0.0);


    // define LSPG problem
    using time_discrete_ops_t = time_discrete_ops;
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_sm_t, ode_case, decoder_t, lspg_state_t, time_discrete_ops_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;

    time_discrete_ops_t tdOps(meshGraph, smGidsToFGidsMap);
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0, tdOps);


    // solvers (linear and GN)
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using hessian_t  = rompp::core::Matrix<eig_dyn_mat>;

    // linear solver
    using solver_tag   = rompp::solvers::linear::iterative::Bicgstab;
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    lin_solver_t linSolverObj;

    // GN solver
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, lin_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-10);
    solver.setMaxIterations(50);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM,
    				t0, dt, Nsteps, solver);

    // // compute the fom corresponding to our rom final state
    // const auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);
  }









  // // do SVD and compute basis
  // Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  // const auto U = svd.matrixU();
  // std::cout << std::setprecision(15) << U << std::endl;
  // // std::cout << S << std::endl;
  // std::cout << "Done with snapshots" << std::endl;

  // // std::ofstream file;
  // // file.open( "phi.txt" );
  // // file << std::fixed
  // //      << std::setprecision(15)
  // //      << U
  // //      << std::endl;
  // // file.close();

  // // read gold basis from file
  // std::vector<std::vector<double>> goldU;
  // readMatrixFromFile("gold_basis.txt", goldU, Nsteps+1);

  // // check that computed matches gold
  // if ( (size_t) goldU.size() != (size_t) U.rows() )
  //   checkStr = "fail";

  // std::cout.precision(15);
  // for (auto i=0; i<U.rows(); i++)
  //   for (auto j=0; j<U.cols(); j++){
  //     auto err = std::abs(goldU[i][j] - U(i,j));
  //     std::cout << "gold = " << goldU[i][j]
  // 		<< " computed = " << U(i,j)
  // 		<< " |err| = " << err
  // 		<< std::endl;

  //     if ( err > eps )
  // 	checkStr = "FAILED";
  //   }

  std::cout << checkStr << std::endl;
  return 0;
}

















// constexpr double eps = 1e-10;
// std::string checkStr {"PASSED"};

// template <typename T>
// void checkSol(const T & y,
// 	      const std::vector<double> & trueS){
//   if (trueS.empty()) {
//     std::cout << " true solution not found, empty " << std::endl;
//     checkStr = "FAILED";
//   }
//   for (size_t i=0; i<trueS.size(); i++){
//     const auto err = std::abs(y[i] - trueS[i]);
//     std::cout << std::fixed << std::setprecision(15)
// 	      << " true = " << trueS[i]
// 	      << " y = " << y[i]
// 	      << " err = " << err
// 	      << std::endl;
//     if ( err > eps or std::isnan(y[i])) checkStr = "FAILED";
//   }
// }

// constexpr bool do_print = true;

// struct Observer{
//   Observer() = default;

//   template <typename T>
//   void operator()(size_t step,
//   		  double t,
//   		  const T & y)
//   {
//     if (do_print){
//       if (step % 50 == 0){
//     	std::ofstream file;
//     	file.open( "sol_" + std::to_string(step) + ".txt" );
//     	for(auto i=0; i < y.size(); i++){
//     	  file << std::fixed << std::setprecision(14) << y[i] ;
//     	  file << std::endl;
//     	}
//     	file.close();
//       }
//     }
//   }
// };


// template <typename T>
// void readMeshGraph(std::string filename,
// 		   std::vector<T> & A){
//   std::ifstream file(filename);
//   std::string line;

//   T lineGIDs;
//   while (getline(file, line)){
//     std::istringstream ss(line);

//     int entry;
//     int j = 0;
//     while (ss >> entry){
//       lineGIDs[j] = entry;
//       j++;
//     }
//     A.emplace_back(lineGIDs);
//   }

//   for (auto & it : A){
//     for (auto & it2 : it){
//       std::cout << " " << it2;
//     }
//     std::cout << std::endl;
//   }
// }

// void readMappingGIDs(std::string filename,
// 		     std::vector<std::array<int,2>> & A){
//   std::ifstream file(filename);
//   std::string line;
//   std::array<int,2> lineGIDs;
//   while (getline(file, line)){
//     std::istringstream ss(line);

//     int entry;
//     int j = 0;
//     while (ss >> entry){
//       lineGIDs[j] = entry;
//       j++;
//     }
//     A.emplace_back(lineGIDs);
//   }

//   for (auto & it : A){
//     for (auto & it2 : it){
//       std::cout << " " << it2;
//     }
//     std::cout << std::endl;
//   }
// }



// struct updateOps{
//   using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
//   using scalar_t	= typename app_t::scalar_type;
//   using s_t		= typename app_t::state_type;
//   using r_t		= typename app_t::residual_type;

//   static void do_update(s_t & v, const scalar_t c,
//   			const r_t & v0, const scalar_t a,
//   			const r_t & v1, const scalar_t b){
//     for (auto i=0; i<v.size(); ++i)
//       v[i] = c*v[i] + a*v0[i] + b*v1[i];
//   }

//   static void do_update(s_t & v,
//   			const r_t & v0, const scalar_t a,
//   			const r_t & v1, const scalar_t b){
//     for (auto i=0; i<v.size(); ++i)
//       v[i] = a*v0[i] + b*v1[i];
//   }

//   static void do_update(s_t & v,
// 			const r_t & v1, const scalar_t & b,
// 			const r_t & v2, const scalar_t & c,
// 			const r_t & v3, const scalar_t & d,
// 			const r_t & v4, const scalar_t & e) {
//     for (auto i=0; i<v.size(); ++i)
//       v[i] = b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
//   }

//   static void do_update(s_t & v, const scalar_t & a,
// 			const r_t & v1, const scalar_t & b,
// 			const r_t & v2, const scalar_t & c,
// 			const r_t & v3, const scalar_t & d,
// 			const r_t & v4, const scalar_t & e) {
//     for (auto i=0; i<v.size(); ++i)
//       v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
//   }
// };

// struct myops{
//   using update_op = updateOps;
// };


// int main(int argc, char *argv[]){
//   using app_t		= rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
//   using scalar_t	= typename app_t::scalar_type;
//   using app_state_t	= typename app_t::state_type;
//   using app_residual_t	= typename app_t::residual_type;
//   using app_jacob_t	= typename app_t::jacobian_type;
//   constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

//   /*
//     mesh.dat: contains list of gid and for each gid the neighbors needed
// 	to compute the residual at that location.
// 	* the neighbors are ordered as east, north, west, south.
// 	* if a neighbor is = -1, it means the current grid point is near boundary
// 	and since neighbors are ordered, we know which boundary is near this point.
// 	* if a grid point is meant to only hold state, then all neighbors are = -2

//     gappy_to_full_gid_map.dat: it is a 2-column file where the first column
// 	contains the gids of each grid point enumerated according to the sample mesh
// 	while the second column contains the corresponding gid of each grid within
// 	the full mesh.
//    */

//   typename app_t::graph_t meshGraph;
//   readMeshGraph("mesh.dat", meshGraph);

//   typename app_t::gids_map_t gToFgidMap;
//   readMappingGIDs("gappy_to_full_gid_map.dat", gToFgidMap);

//   constexpr int Nx = 4, Ny = 4;
//   app_t appobj(Nx, Ny, meshGraph, gToFgidMap);
//   appobj.setup();
//   const auto y0n = appobj.getInitialState();
//   const auto r0n = appobj.residual(y0n, zero);
//   const auto j0nSP = appobj.jacobian(y0n, zero);
//   Eigen::MatrixXd j0n(j0nSP);

//   for(auto i=0; i < r0n.size(); i++)
//     std::cout << i << " "
// 	      << std::fixed << std::setprecision(15)
// 	      << r0n[i] << std::endl;

//   std::cout << "\n";

//   for(auto i=0; i < j0n.rows(); i++){
//     for(auto j=0; j < j0n.cols(); j++){
//       if ( std::abs(j0n(i,j)) < 1e-12)
//   	std::cout << std::fixed << std::setprecision(1) << j0n(i,j) << " ";
//       else
// 	std::cout << std::fixed << std::setprecision(5) << j0n(i,j) << " ";
//     }
//     std::cout << std::endl;
//   }


//   // // if (do_print){
//   //   auto X = appobj.getX(); auto Y = appobj.getY();
//   //   std::ofstream file; file.open( "xy.txt" );
//   //   for(auto i=0; i < X.size(); i++){
//   //     std::cout << std::fixed
//   // 		<< std::setprecision(15)
//   // 		<< X[i] << " " << Y[i] << " "
//   // 		<< y0n[i] << " "
//   // 		<< r0n[i*4] << " " << r0n[i*4+1] << " "
//   // 		<< r0n[i*4+2] << " " << r0n[i*4+3];
//   //     std::cout << std::endl;

//   //     file << std::setprecision(14)
//   // 	   << X[i] << " " << Y[i];
//   //     file << std::endl;
//   //   }
//   //   file.close();
//   // // }

//   // ////////////////////////////////////////
//   // using ode_state_t = rompp::core::Vector<app_state_t>;
//   // using ode_res_t   = rompp::core::Vector<app_residual_t>;
//   // ode_state_t y(y0n);
//   // ode_res_t r(r0n);

//   // constexpr auto ode_case = rompp::ode::ExplicitEnum::RungeKutta4;
//   // using stepper_t = rompp::ode::ExplicitStepper<
//   //   ode_case, ode_state_t, app_t, ode_res_t, scalar_t, myops>;
//   // stepper_t stepperObj(y, appobj, r);

//   // // integrate in time
//   // constexpr scalar_t dt = 0.00001;
//   // constexpr auto Nsteps = 1000;
//   // constexpr scalar_t fint = Nsteps*dt;
//   // Observer obs;
//   // rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs);
//   // //std::cout << std::fixed << std::setprecision(15) << *y.data() << std::endl;

//   // for(auto i=0; i < X.size(); i++){
//   //   std::cout << std::fixed
//   // 	      << std::setprecision(15)
//   // 	      << X[i] << " " << Y[i] << " "
//   // 	      << (*y.data())[i*4] << " "
//   // 	      << (*y.data())[i*4+1] << " "
//   // 	      << (*y.data())[i*4+2] << " "
//   // 	      << (*y.data())[i*4+3];
//   //   std::cout << std::endl;
//   // }
//   ////////////////////////////////////////


//   // using ode_state_t = rompp::core::Vector<app_state_t>;
//   // using ode_res_t   = rompp::core::Vector<app_residual_t>;
//   // using ode_jac_t   = rompp::core::Matrix<app_jacob_t>;

//   // ode_state_t y(y0n);
//   // constexpr auto ode_case = rompp::ode::ImplicitEnum::Euler;
//   // using stepper_t = rompp::ode::ImplicitStepper<
//   //   ode_case, ode_state_t, ode_res_t, ode_jac_t, app_t>;
//   // stepper_t stepperObj(y, appobj);

//   // // define solver
//   // using lin_solver_t = rompp::solvers::iterative::EigenIterative<
//   //   rompp::solvers::linear::iterative::Bicgstab, ode_jac_t>;
//   // rompp::solvers::NewtonRaphson<scalar_t, lin_solver_t> solverO;
//   // solverO.setTolerance(1e-6);
//   // solverO.setMaxIterations(200);

//   // // integrate in time
//   // constexpr scalar_t dt = 0.0001;
//   // constexpr auto Nsteps = 10;
//   // constexpr scalar_t fint = Nsteps*dt;
//   // Observer obs;
//   // rompp::ode::integrateNSteps(stepperObj, y, 0.0, dt, Nsteps, obs, solverO);
//   // std::cout << std::fixed << std::setprecision(14) << *y.data() << std::endl;
//   // {
//   //   using namespace rompp::apps::test;
//   //   checkSol(y,
//   // 	     NonLinAdvDiffReacFlame2dImpGoldStates<ode_case>::get(Nx, Ny, dt, fint));
//   // }


//   // {
//   //   using namespace rompp::apps::test;
//   //   checkSol(y,
//   // 	     NonLinAdvDiffReacFlame2dExpGoldStates<ode_case>::get(Nx,
//   // 								  Ny,
//   // 								  dt,
//   // 								  fint));
//   // }

//   std::cout << checkStr << std::endl;
//   return 0;
// }
