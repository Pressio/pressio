template <typename system_type, typename state_type, typename hessian_type, typename gradient_type, typename linear_solver_type>
class my_gauss_newton_class{
public:
  using scalar_type = typename system_type::scalar_type;
  hessian_type hessian_ = {};
  gradient_type gradient_ = {};
  gradient_type dx_ = {};
  linear_solver_type & linear_solver_ = {};

  my_gauss_newton_class(const system_type & system,const state_type & stateIn,linear_solver_type & linear_solver) : hessian_(system.createHessianObject(stateIn)) , gradient_(system.createGradientObject(stateIn)),dx_(system.createGradientObject(stateIn)), linear_solver_(linear_solver){}; 

  void my_gauss_newton(const system_type & sys, state_type & state, int romSize, int numStepsInWindow)
  {

    double gnorm = 1.;
    double tol = 1e-8;
    scalar_type rnorm = 0.; 
    constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    int iteration = 0;
    //const int sz = romSize*numStepsInWindow;
//    scalar_type * data = eigMatrix.data();
//    Eigen::VectorXd eigVector(data);
//    ::pressio::containers::Vector<Eigen::VectorXd> pressioVector(eigVector);
//    eigenVector.data() = data;  
    ::pressio::containers::Vector<Eigen::Map<Eigen::Matrix<scalar_type, -1, 1,Eigen::ColMajor>>> stateV((*state.data()).data(), (*state.data()).size());

    while (gnorm >= tol){
      sys.computeHessianAndGradient(state,hessian_,gradient_); 
      linear_solver_.solve(hessian_,gradient_,dx_);
      iteration += 1;
      gnorm = (*gradient_.data()).squaredNorm();
      std::cout << "grad norm = " << gnorm << std::endl; 
      std::cout << "Hess norm = " << (*hessian_.data()).squaredNorm() << std::endl; 
      std::cout << "x    norm = " << (*dx_.data()).squaredNorm() << std::endl; 
      std::cout << "iteration = " << iteration << std::endl; 
      hessian_.setZero();
      (*stateV.data()) = (*stateV.data()) + (*dx_.data());
    } 

 } 
};

