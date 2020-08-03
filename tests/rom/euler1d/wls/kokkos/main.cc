
#include "pressio_apps.hpp"
#include "pressio_rom.hpp"
//namespace pressio{ namespace apps{ namespace euler1d{ 


template< typename matrix_t, typename fom_t>
class advanceSolRom{
  Kokkos::View<double*> xhat0_;
  Kokkos::View<double*> velocity_hat_;
  Kokkos::View<double*> velocity_;
  Kokkos::View<double*> F_;
  Kokkos::View<double*> u_;
  matrix_t & Phi_;

  int N_cell_;
  int romSize_;
  double dx_;
  fom_t appObj_;
public:

  advanceSolRom(int romSize, matrix_t &Phi, int N_cell, double dx,const fom_t & appObj) 
    : romSize_(romSize),Phi_(Phi), xhat0_("xhat0",romSize), u_("u",N_cell*3), velocity_hat_("v",romSize), N_cell_(N_cell), dx_(dx),
      velocity_("v",N_cell*3), F_("F",(N_cell+1)*3), appObj_(appObj){}
  void operator()( Kokkos::View<double*>  & xhat, double dt){
    double rk4const[ 4 ];
    rk4const[0] = 1./4.;
    rk4const[1] = 2./4.;
    rk4const[2] = 3./4.;
    rk4const[3] = 1.;
    Kokkos::deep_copy( xhat0_, xhat );


    for (int k = 0; k < 4; k++){
      KokkosBlas::gemv("N",1.,Phi_,xhat,0.,u_);
      //pressio::apps::euler1d::computeVelocity(velocity_,F_,u_,N_cell_,dx_);
      appObj_.velocity(u_,0.,velocity_);
      KokkosBlas::gemv("T",1.,Phi_,velocity_,0.,velocity_hat_);
      //xhat = xhat0_ + dt*rk4const[k]*velocity_hat_;
      for (int i=0; i < romSize_;  i++){
        xhat(i) = xhat0_(i) + dt*rk4const[k]*velocity_hat_(i);
      }
    }
  }

};





int main( int argc, char* argv[] )
{

  int N_cell = 500;         // number of rows 2^12
  double L = 1.;
  double dx = L/(N_cell);
  constexpr double dt = 1e-3;
  double et = 1.;
  double t = 0;
  // Check sizes.
  Kokkos::initialize( argc, argv );
  {
  typedef Kokkos::LayoutLeft   Layout;
  typedef Kokkos::View<double*>   ViewVectorType;
  typedef Kokkos::View<double**,Layout>  ViewMatrixType;

  ViewVectorType u0( "u0", 3*N_cell);
  ViewVectorType u( "u", 3*N_cell);
  ViewVectorType F( "F", 3*(N_cell+1));
  ViewVectorType velocity( "v", 3*N_cell);
  ViewVectorType x( "x", N_cell + 1 );
  ViewVectorType::HostMirror h_u = Kokkos::create_mirror_view( u );
  ViewVectorType::HostMirror h_x = Kokkos::create_mirror_view( x );

  // Initialize y vector on host.
  for ( int i = 0; i < N_cell; ++i ) {
    double gamma = 1.4;

    double rhoL = 1;
    double pL = 1.; 

    double rhoR = 0.125;
    double pR = 0.1;

    if (i < static_cast<int>(N_cell/2) ){
      h_u(::pressio::apps::euler1d::index_map(0,i))  = rhoL;
      h_u(::pressio::apps::euler1d::index_map(1,i))  = 0.;
      h_u(::pressio::apps::euler1d::index_map(2,i)) = pL/(gamma - 1.) ; // + 0.5*rhoU^2/rho 
    }
    else{
      h_u(::pressio::apps::euler1d::index_map(0,i))  = rhoR;
      h_u(::pressio::apps::euler1d::index_map(1,i))  = 0.;
      h_u(::pressio::apps::euler1d::index_map(2,i)) = pR/(gamma - 1.) ; // + 0.5*rhoU^2/rho 
    }
  }

  // Initialize x vector on host.
  h_x(0) = 0.;
  for ( int i = 1; i <= N_cell; ++i ) {
    h_x( i ) = h_x(i-1) + dx;
  }

  // Deep copy host views to device views.
  Kokkos::deep_copy( u, h_u );
  Kokkos::deep_copy( x, h_x );

  // Timer products.
  Kokkos::Timer timer;
  double rk4const[ 4 ];
  rk4const[0] = 1./4.;
  rk4const[1] = 2./4.;
  rk4const[2] = 3./4.;
  rk4const[3] = 1.;

  int save_freq = 8;
  int counter = 0;
  /*
  std::ofstream myfile ("solution.txt"); 
  while (t <= et - dt/2.){
    if (counter%save_freq == 0){
      for (int i = 0;i < N_cell; i++){
        myfile << h_u(::pressio::apps::euler1d::index_map(0,i)) << std::endl;
        myfile << h_u(::pressio::apps::euler1d::index_map(1,i)) << std::endl;
        myfile << h_u(::pressio::apps::euler1d::index_map(2,i)) << std::endl;
      }
    }
    Kokkos::deep_copy( u0, u );
  
    for (int k = 0; k < 4; k++){
      pressio::apps::euler1d::computeVelocity(velocity,F,u,N_cell,dx); 
      for (int i=0; i < N_cell; i++){for (int j = 0; j < 3 ; j++){
          u(::pressio::apps::euler1d::index_map(j,i)) = u0(::pressio::apps::euler1d::index_map(j,i)) + dt*rk4const[k]*velocity(::pressio::apps::euler1d::index_map(j,i));
      }}
    }
    t += dt;
    counter += 1;
    std::cout << t << std::endl;
  }
  myfile.close(); 
  double time = timer.seconds();
  printf( "  time( %g s )\n", time);
  */

  int romSize = 61;
  double tmp_val;
  std::ifstream file("basis.dat");
  ViewMatrixType Phi("Phi",N_cell*3,romSize);
  for (int i =0; i < N_cell ; i++){
    for (int j=0 ; j < 3 ; j++){
      for (int k=0; k < romSize ; k++){
       file >> tmp_val;
       //std::cout << tmp_val << std::endl;
       Phi(::pressio::apps::euler1d::index_map(j,i),k)  =tmp_val;}}}

  file.close(); 


  using fom_t = pressio::apps::euler1d::PressioInterface<ViewVectorType,ViewVectorType,ViewVectorType,double>;
  using scalar_t        = typename fom_t::scalar_type;
  using native_state_t  = typename fom_t::state_type;
  using native_dmat_t      = typename fom_t::dense_matrix_type;
  using fom_state_t     = pressio::containers::Vector<native_state_t>;
  using decoder_jac_t	= pressio::containers::MultiVector<native_dmat_t>;
  using rom_state_t	= pressio::containers::Vector<native_state_t>;
  using hessian_t	= pressio::containers::Matrix<native_dmat_t>; 
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
  using ode_tag   = ::pressio::ode::implicitmethods::Euler;

  fom_t appObj(N_cell,dx,romSize);

  const fom_state_t fomStateInitCond(u);
  const fom_state_t fomStateReference(u);
  pressio::ops::set_zero(fomStateReference);
  pressio::rom::wls::window_size_t numStepsInWindow = 10;
  pressio::rom::wls::rom_size_t wlsSize = romSize*numStepsInWindow;
  scalar_t finalTime = 1.0;
  pressio::rom::wls::window_size_t numWindows = (finalTime/dt)/numStepsInWindow;
  decoder_t decoderObj(Phi);// = readBasis<decoder_t>(appObj,ode_tag(),romSize,fomSize);

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
  ViewVectorType xhat( "xhat", romSize);
  const auto & decoderJac = decoderObj.getReferenceToJacobian();

  rom_state_t wlsStateIC(romSize);
  ::pressio::rom::utils::set_gen_coordinates_L2_projection<scalar_t>(linear_solver, decoderJac,
                      fomStateInitCond, fomStateReference,wlsStateIC);

  // NL solver
  /*
  using gn_t = pressio::solvers::nonlinear::composeLevenbergMarquardt_t<
      wls_system_t, pressio::solvers::nonlinear::LMDefaultUpdate,
      pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
      linear_solver_t>;
  */
  using gn_t = pressio::solvers::nonlinear::composeGaussNewton_t<
    wls_system_t,
    pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::ConvergedWhenRelativeGradientNormBelowTol,
    linear_solver_t>;

  double res;  
  res = KokkosBlas::dot( *wlsState.data(), *wlsState.data() );
  //::pressio::ops::dot(wlsState,wlsState);
  gn_t GNSolver(wlsSystem, wlsState, linear_solver);
  GNSolver.setTolerance(1e-8);
  GNSolver.setMaxIterations(100);


    // set initial guess over window (needed to avoid bad initial guesses that yield NaN)
  for (int i =0; i < numStepsInWindow; i++){
    auto wlsViewAssign = ::pressio::containers::span(wlsState,i*romSize,romSize);
    //auto wlsViewCopy = ::pressio::containers::span(wlsStateIC_,(timeStencilSize_ - 1)*romSize_,romSize_);
    ::pressio::ops::deep_copy(wlsViewAssign, wlsStateIC);
  }


  //  solve wls problem
  auto startTime = std::chrono::high_resolution_clock::now();
  std::ofstream myfile ("wls_rom_solution.txt"); 
  save_freq = 1;
  for (auto iWind = 0; iWind < numWindows; iWind++){
    wlsSystem.advanceOneWindow(wlsState, GNSolver, iWind, dt);
    if (counter%save_freq == 0){
      auto wlsView = ::pressio::containers::span(wlsState,(numStepsInWindow-1)*romSize,romSize);
      KokkosBlas::gemv("N",1.,Phi,*wlsView.data(),0.,u);
      for (int i = 0;i < N_cell; i++){
        myfile << u(::pressio::apps::euler1d::index_map(0,i)) << std::endl;
        myfile << u(::pressio::apps::euler1d::index_map(1,i)) << std::endl;
        myfile << u(::pressio::apps::euler1d::index_map(2,i)) << std::endl;
      }
    }
    counter += 1;
  }
  myfile.close(); 
   /*

  rom_state_t yROM("yRom", romSize);

  // define LSPG type
  using ode_tag  = pressio::ode::implicitmethods::Euler;
  using lspg_problem = pressio::rom::lspg::composeDefaultProblem<
    ode_tag, fom_t, rom_state_t, decoder_t>::type;
  using lspg_stepper_t = typename lspg_problem::lspg_stepper_t;
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
  constexpr auto t0 = zero;
  const auto yRef = *fomStateReference.data();
  lspg_problem lspgProblem(appObj, yRef, decoderObj, yROM,t0);

  // linear solver
  using solver_tag2   = pressio::solvers::linear::direct::getrs;
  using linear_solver_t2 = pressio::solvers::linear::Solver<solver_tag2, hessian_t>;
 linear_solver_t2 linSolverObj;

  // GaussNewton solver
  using nls_t = pressio::solvers::nonlinear::composeGaussNewton_t<
    lspg_stepper_t,
    pressio::solvers::nonlinear::DefaultUpdate,
    pressio::solvers::nonlinear::StopWhenCorrectionNormBelowTol,
    linear_solver_t2>;

 nls_t solver(lspgProblem.getStepperRef(), yROM, linSolverObj);
  solver.setTolerance(1e-14);

  // integrate in time
  pressio::ode::advanceNSteps(lspgProblem.getStepperRef(), yROM, 0.0, dt, 10, solver);

  */    

  }
  Kokkos::finalize();


  return 0;
}
//}}}

