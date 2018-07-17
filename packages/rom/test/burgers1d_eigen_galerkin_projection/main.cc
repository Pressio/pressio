
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_EXP"
#include "apps_burgers1d_eigen.hpp"
#include "experimental/rom_galerkin_implicit_residual_policy.hpp"
#include "experimental/rom_galerkin_implicit_jacobian_policy.hpp"
#include "observer.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>


struct mysizer{
 using state_t = core::vector<apps::burgers1dEigen::state_type>;
 static size_t getSize(state_t & obj){
   return obj.size();
 };
  static void matchSize(const state_t & src, state_t & obj){
    obj.resize(src.size());
 };
};

template<typename T>
void printSol(std::string mess, const T & y){
  std::cout << mess << std::endl;
  for (int i=0; i<y.size(); ++i)
    std::cout << std::setprecision(14) << y[i]  << " ";
  std::cout << std::endl;
}


Eigen::MatrixXd readPhi(int nr, int nc)
{
  Eigen::MatrixXd phi;
  phi.resize(nr,nc);
  
  std::ifstream source;
  source.open("bas.txt", std::ios_base::in);
  std::string line;
  int row = 0;
  while (std::getline(source, line) ){
    //make a stream for the line itself
    std::istringstream in(line);
    // tmp variable to store each entry of the file
    std::vector<std::string> cols(nc);
    for (int i=0; i<nc; i++){
      in >> cols[i];
      phi(row, i) = atof(cols[i].c_str());
    }
    row++;
  }
  source.close();
  return phi;
}//end 



int main(int argc, char *argv[])
{
  using native_state_t = apps::burgers1dEigen::state_type;
  using native_jac_t = apps::burgers1dEigen::jacobian_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;

  int numCell = 100;
  
  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, numCell);
  appObj.setup();

  // wrap with core structures
  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using residual_t = state_t;
  native_state_t y0n = appObj.getInitialState();
  // init state vector
  state_t y0(y0n);

  scalar_t dt = 0.01;
  scalar_t final_t = dt*500;
  auto numSteps = static_cast<unsigned int>(final_t/dt);

  //-----------------------------------------------
  // SOLVERS
  //-----------------------------------------------
  // linear solver
  using lin_solve_t =
    solvers::experimental::linearSolver<jac_t, state_t, state_t>;
  lin_solve_t ls;
  // nonlinear solver
  using nonlin_solve_t =
    solvers::experimental::newtonRaphson<state_t, state_t,
  					 jac_t, lin_solve_t>;
  nonlin_solve_t nonls(ls);

  // //-------------------------------
  // // collect snapshots usinf FOM
  // //-------------------------------
  // state_t y(y0);
  // snapshot_collector collObj(numCell, numSteps);

  // using stepper_t = ode::implicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer, nonlin_solve_t>;
  // stepper_t stepperObj(appObj, nonls);
  
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, numSteps, collObj);
  // printSol("", y);
  // // //  collObj.printAll();

  // //-------------------------------
  // // SVD
  // //-------------------------------
  // Eigen::JacobiSVD<Eigen::MatrixXd> svd(collObj.snapshots_, Eigen::ComputeThinU);
  // auto phi_nat = svd.matrixU();
  // auto phiT_nat = phi_nat;
  // phiT_nat.transposeInPlace();
  // std::cout << "phiSize " << phi_nat.rows() << " " << phi_nat.cols() << std::endl;
  // // std::cout << phi_nat << std::endl;
  // // std::cout  << std::endl;
  // // std::cout << phiT_nat << std::endl;

  // Eigen::MatrixXd phi_nat = readPhi(numCell, 10);
  // auto phiT_nat = phi_nat;
  // phiT_nat.transposeInPlace();
  // //  std::cout << std::setprecision(14) << phi_nat << std::endl;
  
  // //-------------------------------
  // // Galerkin ROM
  // //-------------------------------
  // {
  //   core::matrix<decltype(phi_nat)> phi(phi_nat);
  //   core::matrix<decltype(phiT_nat)> phiT(phiT_nat);

  //   // project initial condition
  //   auto y0nr = *phiT.data() * (*y0.data());
  //   state_t y0r(y0nr);
  //   state_t y2(y0nr);
  //   std::cout << "y0r_Size " << y2.size() << std::endl;
  //   // set to zero because we are integrating the increment wrt y0
  //   //y2.setZero();

  //   using res_pol_t = rom::exp::romGalerkinImplicitResidualPolicy<
  //     state_t, residual_t, model_eval_t, scalar_t, mysizer, decltype(phi)>;
  //   res_pol_t resObj(y0r, phi, phiT);
  //   using jac_pol_t = rom::exp::romGalerkinImplicitJacobianPolicy<
  //     state_t, jac_t, model_eval_t, scalar_t, mysizer, decltype(phi)>;
  //   jac_pol_t jaObj(y0r, phi, phiT);
    
  //   using stepper_t = ode::implicitEulerStepper<
  //     state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //     scalar_t, mysizer, nonlin_solve_t, res_pol_t, jac_pol_t>;
  //   stepper_t stepperObj(appObj, nonls, resObj, jaObj);

  //   ode::integrateNSteps(stepperObj, y2, 0.0, dt, numSteps);
  //   auto gg = y2;// + y0r;
  //   state_t tmp( *phi.data() * (*gg.data()) ) ;
  //   printSol("", tmp);
  // }
  // //-------------------------------

  
  // using stepper_t = ode::implicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer, nonlin_solve_t, res_pol_t, jac_pol_t>;
  // stepper_t stepperObj(appObj, nonls, resObj, jaObj);
  
  // snapshot_collector collObj;
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, numSteps, collObj);
  // printSol("", y+y0);

  // //using stepper_t = ode::explicitRungeKutta4Stepper<
  // using stepper_t = ode::explicitEulerStepper<
  //   state_t, residual_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer>;//, res_pol_t>;
  // stepper_t stepperObj(appObj);//, resObj);

  // // // integration details
  // scalar_t dt = 0.01;
  // scalar_t final_t = 35.;//dt*100;
  // snapshot_collector collObj;
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, final_t/dt, collObj);    
  // printSol("", y);
       
  return 0;
}
