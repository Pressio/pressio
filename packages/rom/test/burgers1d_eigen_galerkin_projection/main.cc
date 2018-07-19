
#include <iostream>
#include <iomanip>
#include <fstream>
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_EXP"
#include "apps_burgers1d_eigen.hpp"
#include "observer.hpp"
#include "experimental/rom_galerkin_implicit_residual_policy.hpp"
#include "experimental/rom_galerkin_implicit_jacobian_policy.hpp"
#include "experimental/basis_operator/rom_basis_operator_default.hpp"
#include "experimental/sampling_operator/rom_sampling_operator_identity.hpp"
#include "experimental/weighting_operator/rom_weighting_operator_identity.hpp"


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


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  // define native types
  using native_state_t = apps::burgers1dEigen::state_type;
  using native_jac_t = apps::burgers1dEigen::jacobian_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;
  // define our wrapper types
  using state_t = core::vector<native_state_t>;
  using residual_t = state_t;
  using jac_t = core::matrix<native_jac_t>;

  // create app object
  int numCell = 100; // number of fv cells
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, numCell);
  appObj.setup();

  // wrap with core structures
  auto y0n = appObj.getInitialState();
  state_t y0(y0n);

  // //-----------------------------------------------
  // // SOLVERS
  // //-----------------------------------------------
  // using lin_solve_t
  //   = solvers::experimental::linearSolver<jac_t, state_t, state_t>;
  // lin_solve_t ls;
  // using nonlin_solve_t
  //   = solvers::experimental::newtonRaphson<state_t, state_t, jac_t, lin_solve_t>;
  // nonlin_solve_t nonls(ls);

  // //-------------------------------
  // // SVD
  // //-------------------------------
  // using native_basis_type = Eigen::MatrixXd;
  // using basis_type = core::matrix<native_basis_type>;
  // auto phi_nat = readPhi(numCell, 10);
  // basis_type phi(phi_nat);
  // // auto phiT_nat = phi_nat;
  // // phiT_nat.transposeInPlace();
  // // core::matrix<decltype(phi_nat)> phi(phi_nat);
  // // core::matrix<decltype(phiT_nat)> phiT(phiT_nat);
  // // //  std::cout << std::setprecision(14) << phi_nat << std::endl;

  // using phiOp_t = rom::experimental::basisOperatorDefault<basis_type>;
  // phiOp_t phiOp(phi);
  // // using POp_t = rom::experimental::samplingOperatorIdentity<basis_type>;
  // // POp_t P;
  // using scale_op_t = rom::experimental::weightingOperatorIdentity<jac_t,phiOp_t>;
  // scale_op_t WOp(phiOp, phi.rows(), phi.cols());

  // //-------------------------------
  // // Galerkin ROM
  // //-------------------------------
  // size_t redSize = phi.cols();
  // // project initial condition
  // state_t yr(redSize);
  // phiOp.project(y0, yr);
  // // std::cout << "y " << *y.data() << std::endl;
  
  // //auto y0nr = *phiT.data() * (*y0.data());
  // // state_t y0r(y0nr);
  // // state_t y2(y0nr);
  // // std::cout << "y0r_Size " << y2.size() << std::endl;
  // // // set to zero because we are integrating the increment wrt y0
  // // //y2.setZero();

  // using res_pol_t = rom::exp::romGalerkinImplicitResidualPolicy<
  //   state_t, residual_t, model_eval_t, scalar_t,
  //   mysizer, phiOp_t, scale_op_t>;
  // res_pol_t resObj(y0, yr, phiOp, WOp);
  // using jac_pol_t = rom::exp::romGalerkinImplicitJacobianPolicy<
  //   state_t, jac_t, model_eval_t, scalar_t, mysizer, phiOp_t, scale_op_t>;
  // jac_pol_t jaObj(y0, yr, phiOp, WOp);
    
  // using stepper_t = ode::implicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer, nonlin_solve_t, res_pol_t, jac_pol_t>;
  // stepper_t stepperObj(appObj, nonls, resObj, jaObj);

  // scalar_t dt = 0.01;
  // scalar_t final_t = dt*1;
  // auto numSteps = static_cast<unsigned int>(final_t/dt);  
  // ode::integrateNSteps(stepperObj, yr, 0.0, dt, numSteps);
  // state_t yrFin(y0);
  // phiOp.leftMultiply(yr, yrFin);
  // // state_t tmp( *phi.data() * (*gg.data()) ) ;
  // printSol("", yrFin);
  // //-------------------------------

  return 0;
}





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
