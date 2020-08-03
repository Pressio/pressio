
#include "pressio_apps.hpp"
#include "pressio_rom.hpp"


int main( int argc, char* argv[] )
{

  int N_cell = 500;         // number of rows 2^12
  double L = 1.;
  double dx = L/(N_cell);
  double dt = 1e-3;
  double et = 1.;
  double t = 0;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  std::string checkStr {"PASSED"};
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));
    using fom_t = pressio::apps::euler1d::tpetra::PressioInterface<double>;
    using scalar_t        = typename fom_t::scalar_type;
    using native_state_t  = typename fom_t::state_type;
    using native_dmat_t      = typename fom_t::dense_matrix_type;
    using fom_state_t     = pressio::containers::Vector<native_state_t>;
    using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;

    using execution_space = Kokkos::DefaultExecutionSpace;
    using kll		= Kokkos::LayoutLeft;
    using k1dLl_d		= Kokkos::View<scalar_t*, kll, execution_space>;
    using k2dLl_d		= Kokkos::View<scalar_t**, kll, execution_space>;
    using rom_state_t	= pressio::containers::Vector<k1dLl_d>;
    using hessian_t	= pressio::containers::Matrix<k2dLl_d>;
    using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
    using ode_tag   = ::pressio::ode::implicitmethods::Euler;

    int romSize = 61;
    fom_t appObj(N_cell,dx,romSize, Comm);

    fom_t::state_t U = appObj.getShockTubeIC();
    fom_t::state_t V = appObj.createVelocity();
    fom_t::state_t U0 = appObj.createVelocity();
    /*
    double tmp_val;
    std::ifstream file("basis.txt");
    auto gridMap = appObj.getGridMap();
    fom_t::dense_matrix_type Phi(*gridMap,3,romSize);
    for (int i =0; i < N_cell ; i++){
      for (int j=0 ; j < 3 ; j++){
        for (int k=0; k < romSize ; k++){
          double * PhiView;
          Phi.getLocalRowView(i,k,PhiView);
          file >> PhiView[j];}}}
  
    file.close(); 


    const fom_state_t fomStateInitCond(U);
    fom_state_t fomStateReference(U);
    pressio::ops::set_zero(fomStateReference);
    pressio::rom::wls::window_size_t numStepsInWindow = 10;
    pressio::rom::wls::rom_size_t wlsSize = romSize*numStepsInWindow;
    scalar_t finalTime = dt*10*2;
    pressio::rom::wls::window_size_t numWindows = (finalTime/dt)/numStepsInWindow;
    decoder_t decoderObj(Phi);

    // build reconstructor for post processing
    using fom_state_reconstr_t    = ::pressio::rom::FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;
    const fom_state_reconstr_t fomStateReconstructor(fomStateReference,decoderObj);


    //  lin solver
    using lin_solver_tag  = pressio::solvers::linear::direct::potrsL;
    using linear_solver_t = pressio::solvers::linear::Solver<lin_solver_tag, hessian_t>;
    linear_solver_t linear_solver;

    //WLS problem ***
    using precon_type = ::pressio::rom::wls::preconditioners::NoPreconditioner;
    using jacobians_update_tag = ::pressio::rom::wls::FrozenJacobian;
    using policy_t     = pressio::rom::wls::HessianGradientSequentialPolicy<fom_t ,  decoder_t,ode_tag,  pressio::matrixLowerTriangular , precon_type,jacobians_update_tag>;

    using wls_system_t = pressio::rom::wls::SystemHessianAndGradientApi<rom_state_t, decoder_t, ode_tag, hessian_t, policy_t>;

    // create policy and wls system
    int jacobianUpdateFrequency = 5;
    policy_t hgPolicy(romSize, numStepsInWindow, decoderObj, appObj, fomStateReference, wls_system_t::timeStencilSize_,jacobianUpdateFrequency);
    wls_system_t wlsSystem(romSize, numStepsInWindow, decoderObj, hgPolicy, fomStateInitCond, fomStateReference, linear_solver);

    // create the wls state
    rom_state_t  wlsState(wlsSize);

    // NL solver
    // 
    //using gn_t = pressio::solvers::nonlinear::composeLevenbergMarquardt_t<
    //    wls_system_t, pressio::solvers::nonlinear::LMDefaultUpdate,
    //    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    //    linear_solver_t>;
    //
    using gn_t = pressio::solvers::nonlinear::composeGaussNewton_t<
      wls_system_t,
      pressio::solvers::nonlinear::DefaultUpdate,
      pressio::solvers::nonlinear::ConvergedWhenRelativeGradientNormBelowTol,
      linear_solver_t>;

    gn_t GNSolver(wlsSystem, wlsState, linear_solver);
    GNSolver.setTolerance(1e-8);
    GNSolver.setMaxIterations(100);


    //  solve wls problem
    auto startTime = std::chrono::high_resolution_clock::now();
    fom_state_t fomStatePP(U);
    std::ofstream myfile ("wls_rom_solution.txt"); 
    int save_freq = 1;
    int counter = 0;
    for (auto iWind = 0; iWind < numWindows; iWind++){
      wlsSystem.advanceOneWindow(wlsState, GNSolver, iWind, dt);
      if (counter%save_freq == 0){
        auto wlsView = ::pressio::containers::span(wlsState,(numStepsInWindow-1)*romSize,romSize);
        fomStateReconstructor(wlsView, fomStatePP);
        auto fomStatePP_native = *fomStatePP.data();
        for (int i = 0;i < N_cell; i++){
          double * dataView;
          fomStatePP_native.getLocalRowView(i,dataView);
          myfile << dataView[0] << std::endl;
          myfile << dataView[1] << std::endl;
          myfile << dataView[2] << std::endl;
        }
      }
      counter += 1;
  }
  myfile.close(); 
  double wlsNormGold = 1.342267626e+02;
  auto wlsNorm = ::pressio::ops::norm2(wlsState);
  std::cout << "WLS State Norm = " << std::setprecision(9) << wlsNorm << std::endl;
  if (std::abs(wlsNorm - wlsNormGold) >= 1e-8){
    checkStr = "FAILED";  
  }
  std::cout << checkStr << std::endl; 
  */
    // This code can be used to run the FOM with rk4 and save to solution.txt 
    double rk4const[ 4 ];
    rk4const[0] = 1./4.;
    rk4const[1] = 2./4.;
    rk4const[2] = 3./4.;
    rk4const[3] = 1.;
    int save_freq = 8;
    int counter = 0;
    //std::ofstream myfile ("solution.txt"); 
    while (t <= et - dt/2.){
      fom_t::state_t U0(U,Teuchos::DataAccess::Copy);
      if (counter%save_freq == 0){
        for (int i = 0;i < N_cell; i++){
          double *uLocal;
          //U.getLocalRowView(i,uLocal);
          //myfile << uLocal[0] << std::endl;
          //myfile << uLocal[1] << std::endl;
          //myfile << uLocal[2] << std::endl;
        }
      }
    
      for (int k = 0; k < 4; k++){
        appObj.velocity(U,0.,V);
        for (int i=0; i < N_cell; i++){
          double *Up;
          double *U0p;
          double *Vp;
          U.getLocalRowView(i,Up);
          U0.getLocalRowView(i,U0p);
          V.getLocalRowView(i,Vp);
          for (int j=0;j<3;j++){
            Up[j] = U0p[j] + dt*rk4const[k]*Vp[j];
          } 
        }
      }
      t += dt;
      counter += 1;
      std::cout << t << std::endl;
    }
    //myfile.close(); 
   

  


  }
  return 0;
}
//}}}

