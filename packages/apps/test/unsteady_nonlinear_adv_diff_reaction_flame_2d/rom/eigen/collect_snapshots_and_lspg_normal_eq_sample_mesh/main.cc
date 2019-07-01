
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG"
#include "APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D"
#include "../../../fom/gold_states_implicit.hpp"

using scalar_t		= double;
using uint_t		= unsigned int;
using eig_dyn_mat	= Eigen::MatrixXd;
using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
constexpr auto zero	= ::rompp::utils::constants::zero<scalar_t>();
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
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacobian_t	= typename app_t::jacobian_type;
  using ode_state_t = rompp::containers::Vector<app_state_t>;
  using ode_res_t   = rompp::containers::Vector<app_rhs_t>;
  using ode_jac_t   = rompp::containers::Matrix<app_jacobian_t>;

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


// TODO: fix the types, this is not the best way to pass things
struct time_discrete_ops{
  using app_sm_t   = rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;
  using graph_t	   = typename app_sm_t::graph_t;
  using gids_map_t = typename app_sm_t::gids_map_t;
  using state_t	   = typename app_sm_t::state_type;
  using resid_t	   = typename app_sm_t::state_type;
  // the mv_t should be taken from the lspg problem
  using mv_t	   = typename app_sm_t::mv_t;
  const int numSpecies = 4;

  const graph_t & graph_;
  const gids_map_t & gidsMap_;

  time_discrete_ops(const graph_t & graph,
		    const gids_map_t & gidsMap)
    : graph_{graph}, gidsMap_{gidsMap}{
  }

  void time_discrete_euler(resid_t & R,
			   const state_t & yn,
			   const state_t & ynm1,
			   scalar_t dt) const{

    std::cout << "size(R) " << R.rows() << " " << R.cols() << std::endl;
    /*
     * On input R contains the application RHS
     * yn is y_n (i.e. state at time step n)
     * ynm1 is y_n-1 (i.e. state at time step n-1)
     * note that R and y have diffferent sizes!
     * but yn and ynm1 have the same size
     */
    // loop over cells where residual needs to be computed
    for (size_t rPt=0; rPt < graph_.size(); ++rPt)
    {
      // get GID of the current residual cell
      const auto rCellGID  = graph_[rPt][0];

      // distinguish between dofs fro residual points vs state points
      // update the residual using corresponding state at that same cell
      // the dofGID_r
      auto dofGID_r  = numSpecies * rPt;
      // the dofGID_s
      auto dofGID_s  = numSpecies * rCellGID;

      for (auto iDof=0; iDof<numSpecies; iDof++)
	R[dofGID_r+iDof] = yn[dofGID_s+iDof] -
			   ynm1[dofGID_s+iDof] -
			   dt*R[dofGID_r+iDof];
    }
  }

  void time_discrete_jacobian(mv_t & jphi,
			      const mv_t & phi,
			      scalar_t prefactor,
			      scalar_t dt) const
  {
    // loop over cells where residual needs to be computed
    for (size_t rPt=0; rPt < graph_.size(); ++rPt)
    {
      // get GID of the current residual cell
      const auto cellGID  = graph_[rPt][0];

      // distinguish between dofs fro residual points vs state points
      // update the residual using corresponding state at that same cell
      // the dofGID_r
      auto dofGID_r  = numSpecies * rPt;
      auto dofGID_s  = numSpecies * cellGID;

      for (auto iDof=0; iDof<numSpecies; iDof++){
	for (auto j=0; j<jphi.cols(); j++)
	  jphi(dofGID_r+iDof, j) = phi(dofGID_s+iDof, j)
	    - prefactor*dt*jphi(dofGID_r+iDof,j);
      }
    }
  }
};




int main(int argc, char *argv[]){
  std::string checkStr {"PASSED"};

  // set parameters to use for runs
  constexpr int Nx = 10, Ny = 10;
  constexpr scalar_t dt = 0.001;
  constexpr auto Nsteps = static_cast<uint_t>(3);

  //---------------------------------------------------
  // run FOM and get snapshots using full mesh case
  // run FOM and collect snapshots
  FomRunner fom(Nx, Ny);
  const auto fomY = fom.run(dt, Nsteps);
  // get snapshots
  const auto & S = fom.getSnapshots();
  //std::cout << S << std::endl;
  // compute SVD to create basis
  Eigen::JacobiSVD<eig_dyn_mat> svd(S, Eigen::ComputeThinU);
  const auto U = svd.matrixU();
  for (auto i=0; i<U.rows(); ++i){
    if (i % 4 == 0){
      std::cout << std::endl;
      std::cout << i/4 << std::endl;
    }
    std::cout << i << " ";
    for (auto j=0; j<U.cols(); ++j)
      std::cout << std::setprecision(15) << U(i,j) << " ";
    std::cout << std::endl;
  }
  // std::cout << std::setprecision(15) << U << std::endl;
  std::cout << "size(U) " << U.rows() << " " << U.cols() << std::endl;
  std::cout << "Done with snapshots" << std::endl;
  //---------------------------------------------------

  /* run ROM with sample mesh */
  {
    // this is the size of the rom
    constexpr int romSize = Nsteps;

    // app type
    using app_sm_t = rompp::apps::UnsteadyNonLinAdvDiffReacFlame2dSampleMeshEigen;

    // read sample mesh
    typename app_sm_t::graph_t meshGraph;
    readMeshGraph("mesh.dat", meshGraph);
    typename app_sm_t::gids_map_t smGidsToFGidsMap;
    readMappingGIDs("gappy_to_full_gid_map.dat", smGidsToFGidsMap);

    // create app object
    app_sm_t appobj(Nx, Ny, meshGraph, smGidsToFGidsMap);
    appobj.setup();

    // typedefs used for rom
    using lspg_state_t	= rompp::containers::Vector<eig_dyn_vec>;
    using decoder_jac_t	= rompp::containers::MultiVector<eig_dyn_mat>;
    using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

    // use sample mesh to extract only taget rows from basis
    constexpr int numSpecies = 4;
    eig_dyn_mat Usm(smGidsToFGidsMap.size()*numSpecies, romSize);
    // loop over cells holding state
    for (size_t iPt=0; iPt < smGidsToFGidsMap.size(); ++iPt)
    {
      // get GID of the current cell
      const auto cellGID  = smGidsToFGidsMap[iPt][0];
      // get the GID in the full mesh of this cell
      const auto cellGIDinFullMesh = smGidsToFGidsMap[iPt][1];

      auto dofGID_s  = numSpecies * cellGIDinFullMesh;
      auto dofGID  = numSpecies * cellGID;
      std::cout << std::endl;
      std::cout << iPt << " " << dofGID_s << " " << dofGID << " \n";
      for (auto iDof=0; iDof<numSpecies; iDof++)
      {
	std::cout << dofGID+iDof << " ";
    	for (auto j=0; j<U.cols(); j++){
    	  Usm(dofGID+iDof, j) = U(dofGID_s+iDof,j);
    	  std::cout << std::setprecision(15) << Usm(dofGID+iDof, j) << " ";
    	}
    	std::cout << "\n";
      }
    }
    std::cout << "size(Usm) " << " "
	      << Usm.rows() << " "
	      << Usm.cols() << std::endl;

    // create decoder
    decoder_jac_t phi(Usm);
    decoder_t decoderObj(phi);

    // my reference state
    auto yRef = appobj.getInitialState();

    // define ROM state and set to zero
    lspg_state_t yROM(romSize);
    yROM.putScalar(0.0);

    // define object with user-defined ops
    using time_discrete_ops_t = time_discrete_ops;
    time_discrete_ops_t tdOps(meshGraph, smGidsToFGidsMap);

    // define LSPG problem
    using lspg_problem_type = rompp::rom::DefaultLSPGTypeGenerator<
      app_sm_t, ode_case, decoder_t, lspg_state_t, time_discrete_ops_t>;
    using lspg_generator = rompp::rom::LSPGUnsteadyProblemGenerator<lspg_problem_type>;
    lspg_generator lspgProblem(appobj, yRef, decoderObj, yROM, t0, tdOps);

    // solvers (linear and GN)
    // hessian comes up in GN solver, it is (J phi)^T (J phi)
    // rom is solved using eigen, hessian is wrapper of eigen matrix
    using hessian_t  = rompp::containers::Matrix<eig_dyn_mat>;

    // linear solver
    using solver_tag   = rompp::solvers::linear::iterative::Bicgstab;
    using lin_solver_t = rompp::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
    lin_solver_t linSolverObj;

    // GN solver
    using lspg_stepper_t = typename lspg_problem_type::lspg_stepper_t;
    using gnsolver_t   = rompp::solvers::iterative::GaussNewton<
      lspg_stepper_t, lin_solver_t>;
    gnsolver_t solver(lspgProblem.stepperObj_, yROM, linSolverObj);
    solver.setTolerance(1e-6);
    solver.setMaxIterations(30);

    // integrate in time
    rompp::ode::integrateNSteps(lspgProblem.stepperObj_, yROM, t0, dt, 3, solver);

    // compute the fom corresponding to our rom final state
    const auto yFomFinal = lspgProblem.yFomReconstructor_(yROM);

    std::cout << "yROM " << std::endl;
    std::cout << *yROM.data() << std::endl;

    scalar_t rmsErr = zero;
    auto nSM = smGidsToFGidsMap.size();
    for (auto i=0; i<nSM; ++i){
      auto SMgid = smGidsToFGidsMap[i][0];
      auto FMgid = smGidsToFGidsMap[i][1];
      for (auto iDof=0; iDof<4; iDof++){
	const auto err = yFomFinal[SMgid*4+iDof] - fomY[FMgid*4+iDof];
	rmsErr += err*err;
      }
    }
    auto err = std::sqrt(rmsErr/static_cast<scalar_t>(nSM));
    std::cout << "error = " << err << std::endl;
    if (err > 1e-6)
      checkStr = "FAILED";
  }

  std::cout << checkStr << std::endl;
  return 0;
}
