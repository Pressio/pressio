
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>

#include "matrix/core_matrix_eigen.hpp"
#include "vector/core_vector_eigen.hpp"

#include "ode_implicit_euler_stepper.hpp"
#include "ode_euler_stepper.hpp"
#include "ode_rk4_stepper.hpp"
#include "ode_integrate_n_steps.hpp"

#include "apps_burgers1d_eigen.hpp"


struct eigenVectorStateResizer{
  using vec_t = Eigen::VectorXd;
  // this has to be default constructible
  // this will be checked at compile-time by ode package
  void operator()(const vec_t & source, vec_t & dest){
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};



struct snapshot_collector{
  using state_t = apps::burgers1dEigen::state_type;
  using matrix = std::vector<state_t>;  
  matrix snapshots_;
  size_t count_;
  void operator()(size_t step, double t, const state_t & x){
    if (step % 50 ==0 ){
      snapshots_.emplace_back(x);
      count_++;
    }
  }
  size_t getCount(){
    return count_;
  };
  void print(){
    auto const & finalSol = snapshots_.back();
    std::cout << "Final solution " << std::endl;
    for (int i=0; i<finalSol.size(); ++i)
      std::cout << finalSol(i)  << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  } 
  void printToFile(){
    std::ofstream file;
    file.open( "out.txt" );
    for (size_t step=0; step<count_; ++step){
      auto const & sol = snapshots_[step];
      for (int i=0; i< sol.size(); ++i)
	file << std::fixed << std::setprecision(10) << sol(i) << " "; 
      file << std::endl;
    }
    file.close();    
  };
};


struct snapshot_collector_myvec{
  using state_t = core::vector<apps::burgers1dEigen::state_type>;
  using matrix = std::vector<state_t>;  
  matrix snapshots_;
  int count_;
  void operator()(size_t step, double t, const state_t & x){
    if (step % 50 ==0 ){
      snapshots_.emplace_back( *x.view() );
      count_++;
    }
  }

  size_t getCount(){ return count_; };

  void printToFile(){
    std::ofstream file;
    file.open( "outImp.txt" );
    for (size_t step=0; step<count_; ++step){
      auto const & sol = snapshots_[step];
      for (int i=0; i< sol.size(); ++i)
	file << std::fixed << std::setprecision(10) << sol[i] << " "; 
      file << std::endl;
    }
    file.close();    
  };
  
};


struct myVectorStateResizer{
  // this has to be default constructible
  using vec_t = core::vector<Eigen::VectorXd>;
  void operator()(const vec_t & source, vec_t & dest){
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};


int main(int argc, char *argv[])
{
  using app_state_t = apps::burgers1dEigen::state_type;
  using app_jac_t = apps::burgers1dEigen::jacobian_type;

  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  std::cout << mu[0] << " " << mu(1) << std::endl;
  apps::burgers1dEigen appObj(mu, 20);
  appObj.setup();

  app_state_t U = appObj.copyInitialState();
  snapshot_collector snColl;
  ode::eulerStepper<app_state_t,app_state_t,double,eigenVectorStateResizer> myStepper;
  ode::integrateNSteps(myStepper, appObj, U, snColl, 0.0, 0.0035, 10000);
  std::cout << snColl.getCount() << std::endl;
  snColl.print();

  // // wrap state vector
  // using myvec_t = core::vector<app_state_t>;
  // myvec_t y( appObj.copyInitialState() ); // y contains the initial condition of app
  // // wrap jacobian
  // using mymat_t = core::matrix<app_jac_t>;
  // snapshot_collector_myvec coll2;
  // ode::implicitEulerStepper<myvec_t,myvec_t, mymat_t, apps::burgers1dEigen,
  // 			    double, myVectorStateResizer> myStepperImp(appObj);
  // ode::integrateNStepsImpl(myStepperImp, y, coll2, 0.0, 0.07, 500);
  // //std::cout << coll2.getCount() << std::endl;

  return 0;
}
