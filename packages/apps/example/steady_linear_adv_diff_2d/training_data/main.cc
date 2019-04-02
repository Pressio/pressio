
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_STEADY"
#include "APPS_STEADYLINADVDIFF2D"
#include "utils_epetra.hpp"
#include <random>

using gen_t	   = std::mt19937;
using rand_distr_t = std::uniform_real_distribution<double>;

// range of Prandtl number
constexpr std::array<double,2> Pr_range{1.0, 5.0};

// range of Reynolds number
constexpr std::array<double,2> Re_range{10., 100.0};

struct ResidualSampler{
  using vec_t = rompp::core::Vector<Epetra_Vector>;

  mutable std::vector<double> myR_;

  std::vector<double> getR() const{
    return myR_;
  };

  void observeResidualWhenSolverConverged(const vec_t & R) const{
    auto map = R.data()->Map();
    auto N = map.NumMyElements();
    if (myR_.size() != (size_t) N)
      myR_.resize(N);

    for (auto i=0; i<N; i++)
      myR_[i] = R[i];
  }
};


int main(int argc, char *argv[]){
  using true_fom_t	= rompp::apps::SteadyLinAdvDiff2dEpetra;
  using fom_adapter_t	= rompp::apps::SteadyLinAdvDiff2dEpetraRomAdapter;
  using scalar_t	= typename fom_adapter_t::scalar_type;
  using native_state	= typename fom_adapter_t::state_type;

  // MPI
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 1);



  // the discretization to use for solver
  const int Nx = 11, Ny = Nx*2-1;
  // tot num of dof
  const int numDof = (Nx-2)*Ny;

  // TODO hardcode residual sample points for now:
  std::vector<int> residualSampleDofs;

  const int numResSamples = 50;
  residualSampleDofs.resize(numResSamples);
  for (int i; i<numResSamples; i++) residualSampleDofs[i] = i*3;

  // open files for MLEM python scripts

  // training_features.dat
  std::ofstream file_features;
  file_features.open( "training_features.dat" );

  // training_qoi_fom.dat
  std::ofstream file_qoi_fom;
  file_qoi_fom.open( "training_qoi_fom.dat" );

  // training_qoi_rom.dat
  std::ofstream file_qoi_rom;
  file_qoi_rom.open( "training_qoi_rom.dat" );


  //---------------------------------------
  // generate random samples of Parameters

  // random number generator (seeded)
  unsigned int seed = 12;
  std::mt19937 engine(seed);
  rand_distr_t distr(0., 1.0);
  auto genPr = [&distr, &engine](){
		 auto c1 = Pr_range[1]-Pr_range[0];
		 auto c2 = Pr_range[0];
		 return c1 * distr(engine) + c2;
	     };

  auto genRe = [&distr, &engine](){
		 auto c1 = Re_range[1]-Re_range[0];
		 auto c2 = Re_range[0];
		 return c1 * distr(engine) + c2;
	     };

  // number of sample to take
  constexpr int nSamples = 30;

  // fill for Prandtl
  std::vector<double> PrS(nSamples);
  std::generate(PrS.begin(), PrS.end(), genPr);

  // fill for Reynolds
  std::vector<double> ReS(nSamples);
  std::generate(ReS.begin(), ReS.end(), genRe);

  if(rank==0){
    auto it1 = PrS.begin();
    auto it2 = ReS.begin();
    for( ;it1<PrS.end(), it2<ReS.end(); it1++, it2++)
      std::cout << std::setprecision(15)
		<< *it1 << " " << *it2
		<< "\n";
  }





//TODO loop over parameter pairs, solve ROM and FOM, output and save residual
  for (int iSample=0; iSample<nSamples; iSample++)
  {

	  //  {

	  // we run the FOM and LSPG for same values of parameters
	  scalar_t Pr = PrS[iSample];
	  scalar_t Re = ReS[iSample];

	  // -------------------------
	  // -------------------------
	  // run FOM model first
	  // -------------------------
	  // -------------------------

	  true_fom_t  appObj(Comm, Nx, Ny, Pr, Re);
	  appObj.setup();
	  appObj.assembleMatrix();
	  appObj.fillRhs();
	  appObj.solve();
	  //appObj.printStateToFile("fom.txt");
	  rompp::core::Vector<native_state> yFom(*appObj.getState());

	  // -------------------------
	  // -------------------------
	  // do LSPG ROM
	  // -------------------------
	  // -------------------------
	  using native_state	= typename fom_adapter_t::state_type;
	  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
	  using lspg_state_t	= rompp::core::Vector<eig_dyn_vec>;
	  using decoder_jac_t	= rompp::core::MultiVector<Epetra_MultiVector>;
	  using decoder_t	= rompp::rom::LinearDecoder<decoder_jac_t>;

	  constexpr int romSize = 5;

	  // app object for running rom
	  fom_adapter_t  appObjROM(Comm, Nx, Ny, Pr, Re);

	  // store modes from file
	  const decoder_jac_t phi =
			  rompp::apps::test::epetra::readBasis("basis.txt", romSize, numDof,
					  Comm, appObjROM.getDataMap());
	  // decoder object
	  decoder_t decoderObj(phi);

	  // my reference state
	  auto yRef = appObjROM.getState();
	  yRef->PutScalar( static_cast<scalar_t>(0) );

	  // define ROM state and set to zero
	  lspg_state_t yROM(romSize);
	  yROM.putScalar(0.0);

	  // define LSPG type
	  using lspg_problem_type = rompp::rom::DefaultLSPGSteadyTypeGenerator<
			  fom_adapter_t, decoder_t, lspg_state_t>;
	  rompp::rom::LSPGSteadyProblemGenerator<lspg_problem_type> lspgProblem(
			  appObjROM, *yRef, decoderObj, yROM);

	  // create sampler for residual (this is passed to solver)
	  using observer_t	 = ResidualSampler;
	  observer_t myResidSampler;

	  // GaussNewton solver
	  // hessian comes up in GN solver, it is (J phi)^T (J phi)
	  // rom is solved using eigen, hessian is wrapper of eigen matrix
	  using eig_dyn_mat	 = Eigen::Matrix<scalar_t, -1, -1>;
	  using hessian_t	 = rompp::core::Matrix<eig_dyn_mat>;
	  using solver_tag	 = rompp::solvers::linear::LSCG;
	  using converged_when_t = rompp::solvers::iterative::default_convergence;
	  using rom_system_t	 = typename lspg_problem_type::lspg_system_t;
	  using gnsolver_t	 = rompp::solvers::iterative::GaussNewton<
			  scalar_t, solver_tag, rompp::solvers::EigenIterative,
			  converged_when_t, rom_system_t, hessian_t, void, void, void, observer_t>;
	  gnsolver_t solver(lspgProblem.systemObj_, yROM, myResidSampler);
	  solver.setTolerance(1e-14);
	  solver.setMaxIterations(200);
	  solver.solve(lspgProblem.systemObj_, yROM);

      // Write out training data needed for MLEM

	  // Features: parameter values + 50 residual samples

	  // retrieve the residual that was stored by the sampler
	  auto storedR = myResidSampler.getR();
	  // sample randomly


	  // write to file
	  file_features << std::scientific << std::setprecision(14) << Re << " " << Pr << " " ;
	  for(auto i=residualSampleDofs.begin(); i !=residualSampleDofs.end(); i++){
		  file_features << std::scientific << std::setprecision(14) << storedR[*i] << " ";
	  }
	  file_features << std::endl;


	  // ROM quantity of interest
	  auto yFomApprox = lspgProblem.yFomReconstructor_(yROM);
	  const auto energyROM = rompp::core::ops::norm2(yFomApprox);
	  file_qoi_rom << energyROM << std::endl;

	  // FOM quantity of interest
	  const auto energyFOM = rompp::core::ops::norm2(yFom);
      file_qoi_fom << energyFOM << std::endl;


  }
  // Loop Over Parameters

  // close files
  file_features.close();
  file_qoi_fom.close();
  file_qoi_rom.close();


  MPI_Finalize();
  return 0;
}
