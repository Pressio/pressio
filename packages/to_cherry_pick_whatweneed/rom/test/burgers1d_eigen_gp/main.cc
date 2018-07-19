
#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>

#include "ode_euler_stepper.hpp"
#include "ode_implicit_euler_stepper.hpp"
#include "ode_integrate_n_steps.hpp"

#include "apps_burgers1d_eigen.hpp"
#include "rom_gp_draft.hpp"
#include "rom_lspg_draft.hpp"

#include "matrix/core_matrix_traits.hpp"
#include "matrix/core_matrix_eigen.hpp"
#include "vector/core_vector_traits.hpp"
#include "vector/core_vector_serial_arbitrary.hpp"
#include "vector/core_vector_std.hpp"
#include "vector/core_vector_eigen.hpp"

#include "svd_solver_eigen.hpp"
#include "svd_solver_traits.hpp"


struct eigenVectorStateResizer{
  using vec_t = Eigen::VectorXd;
  // this has to be default constructible
  // this will be checked at compile-time by ode package  
  void operator()(const vec_t & source, vec_t & dest)
  {
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};


struct myVectorStateResizer{
  using vec_t = core::Vector<Eigen::VectorXd>;
  // this has to be default constructible
  // this will be checked at compile-time by ode package
  void operator()(const vec_t & source, vec_t & dest)
  {
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};


bool checkTime(double t, double targetT)
{
  return std::abs(t-targetT)<1e-8;
}

struct snapshot_collector{
  using state_t = apps::Burgers1dEigen::state_type;
  using matrix_t = Eigen::MatrixXd;
  matrix_t snapshots_;
  size_t count_;
  int j;
  void operator()(size_t step, double t, const state_t & x)
  {
    if (snapshots_.rows() == 0){
      j = 0;
      snapshots_.resize(x.rows(), 11);
    }   
    if ( checkTime(t, 0.0) || checkTime(t, 3.5) ||
	 checkTime(t, 7.0) || checkTime(t, 10.5) ||
	 checkTime(t, 14.0)|| checkTime(t, 17.5)||
	 checkTime(t, 21.0)|| checkTime(t, 24.5)||
	 checkTime(t, 28.0)|| checkTime(t, 31.5)||
	 checkTime(t, 35.0)
	){
      std::cout << t << std::endl;
      for (int i=0; i<x.rows(); ++i)
	snapshots_(i,j) = x(i);
      count_++;
      j++;
    }
  }

  matrix_t getSnaps() const{ return snapshots_; };
  size_t getCount() const{ return count_; };
  void printe(){
    for (int i=0; i<snapshots_.rows(); ++i){
      for (int k=0; k<snapshots_.cols(); ++k){
	std::cout << snapshots_(i,k) << " ";
      }      
      std::cout << std::endl;
    }
  }
  void printFile(){
    std::ofstream file;
    file.open( "out.txt" );
    for (int i=0; i<snapshots_.rows(); ++i){
      for (int k=0; k<snapshots_.cols(); ++k){
  	file << std::fixed
	     << std::setprecision(10) << snapshots_(i,k) << " ";
      }
      file << std::endl;
    }
    file.close();    
  }

};
//////////////////////////////////

template <typename svdsolve_type>
struct snapshot_collector_myvec{
  using state_t = core::Vector<apps::Burgers1dEigen::state_type>;
  using matrix_t = Eigen::MatrixXd;

  matrix_t snapshots_;
  size_t count_;
  state_t Udoty_;
  using lsv_matrix_type = typename svd::details::traits<svdsolve_type>::u_matrix_t;    
  lsv_matrix_type lsv_;

  snapshot_collector_myvec(svdsolve_type & svd){
    lsv_ = svd.leftSingularVectors();
  }
  
  int j;
  void operator()(size_t step, double t, const state_t & x)
  {
    if (snapshots_.rows() == 0){
      j = 0;
      snapshots_.resize(lsv_.rows(),11);
      Udoty_.resize(lsv_.rows());
    }   
    if ( checkTime(t, 0.0) ||
	 checkTime(t, 3.5) ||
	 checkTime(t, 7.0) ||
	 checkTime(t, 10.5) ||
	 checkTime(t, 14.0)||
	 checkTime(t, 17.5)||
	 checkTime(t, 21.0)||
	 checkTime(t, 24.5)||
	 checkTime(t, 28.0)||
	 checkTime(t, 31.5)||
	 checkTime(t, 35.0)
	){
      std::cout << t << std::endl;
      x.matMultiply(lsv_, Udoty_);
      for (int i=0; i<Udoty_.size(); ++i)
	snapshots_(i,j) = Udoty_[i];
      count_++;
      j++;
    }
  }

  size_t getCount() const{ return count_; };

  void printFile(std::string filename){
    std::ofstream file;
    file.open(filename);
    for (int i=0; i<snapshots_.rows(); ++i){
      for (int k=0; k<snapshots_.cols(); ++k){
  	file << std::fixed
	     << std::setprecision(10) << snapshots_(i,k) << " ";
	std::cout << std::setprecision(10) << snapshots_(i,k) << " ";
      }
      file << std::endl;
      std::cout << std::endl;
    }
    file.close();    
  }
};




int main(int argc, char *argv[])
{
  using app_state_t = apps::Burgers1dEigen::state_type;
  
  //-----------------
  // create the app 
  //-----------------
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  apps::Burgers1dEigen appObj(mu, 200);
  appObj.setup();
  app_state_t U = appObj.copyInitialState();

  //-------------------
  // collect snapshots 
  //-------------------
  auto snColl = std::make_shared<snapshot_collector>();
  ode::eulerStepper<app_state_t,app_state_t,double,eigenVectorStateResizer> myStepper;
  ode::integrateNSteps(myStepper, appObj, U, *snColl, 0.0, 0.07, 501);
  std::cout << snColl->getCount() << std::endl;
  snColl->printFile();
    
  // wrap snaps into matrix (this can be done inside the collector directly)
  using mymat_t = core::Matrix<Eigen::MatrixXd>;
  auto snaps = snColl->getSnaps();
  mymat_t MM(snaps);

  //-------------------
  // do SVD
  //-------------------
  using svd_type = svd::solver<mymat_t, svd::svdKind::EigenJacobi>;
  svd_type svdSolve;
  svdSolve.compute(MM);
  
  //--------------------------------------
  // wrap the app state with our vector
  //--------------------------------------
  using myvec_t = core::Vector<app_state_t>;
  myvec_t y( appObj.copyInitialState() ); // y contains the initial condition of app
  
  // //---------------------
  // // galerkin projection
  // //---------------------
  rom::GP<myvec_t,apps::Burgers1dEigen,svd_type> gpSolver(y, appObj, svdSolve);
  snapshot_collector_myvec<svd_type> snColl2(svdSolve);
  ode::eulerStepper<myvec_t,myvec_t,double,myVectorStateResizer> myStepper2;
  ode::integrateNSteps(myStepper2, gpSolver, y, snColl2, 0.0, 0.07, 501); //0, to 35
  snColl2.printFile("outGP.txt");

  //---------------------
  // lspg
  //---------------------
  myvec_t y2( appObj.copyInitialState() ); // y contains the initial condition of app
  
  using app_jac_t = apps::Burgers1dEigen::jacobian_type;
  using mymatjac_t = core::Matrix<app_jac_t>;

  using rom_algo_type = rom::lspg<myvec_t, mymatjac_t,apps::Burgers1dEigen,svd_type>;
  rom_algo_type lspgObj(y2, appObj, svdSolve);

  ode::ImplicitEulerStepper<myvec_t,myvec_t, mymatjac_t, rom_algo_type,
  			    double, myVectorStateResizer> myStepperImp(lspgObj, true);
  snapshot_collector_myvec<svd_type> snCollEI(svdSolve);
  ode::integrateNStepsImpl(myStepperImp, y2, snCollEI, 0.0, 0.07, 501);
  snCollEI.printFile("outLSPG.txt");
  
  return 0;
}
