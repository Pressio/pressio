
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>

#include "ode_euler_stepper.hpp"
#include "ode_rk4_stepper.hpp"
#include "ode_integrate_n_steps.hpp"
// // #include "vector/core_vector_traits.hpp"
// // #include "vector/core_vector_epetra.hpp"
// // #include "vector/core_vector_serial_arbitrary.hpp"
// // #include "vector/core_vector_std.hpp"
// // #include "vector/core_vector_eigen.hpp"

#include "apps_burgers1d_eigen.hpp"


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
    for (int i=0; i<finalSol.size(); ++i)
      std::cout << finalSol(i)  << std::endl;
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



int main(int argc, char *argv[])
{
  using state_t = apps::burgers1dEigen::state_type;

  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  std::cout << mu[0] << " " << mu(1) << std::endl;
  apps::burgers1dEigen appObj(mu);
  appObj.setup();

  state_t U = appObj.copyInitialState();
  snapshot_collector snColl;
  ode::eulerStepper<state_t,state_t,double,eigenVectorStateResizer> myStepper;
  ode::integrateNSteps(myStepper, appObj, U, snColl, 0.0, 0.0035, 10000);
  std::cout << snColl.getCount() << std::endl;
  snColl.print();
  
  return 0;
}





// template <typename state_type, typename oapp>
// class GP()
// {
// private:
//   state_type myY_;
//   oapp * appPtr_;
    
// public:
//   GP(state_type & src, ...)
//     : myY(src), ... {}
//   ~GP(){}
  
//   void operator()()
//   {
//     (*oappPtr)()( V' * y )
//   }

//   void run()
//   {
//     ode::eulerStepper<state_t,state_t,double,stateResizer> myStepper;
//     ode::integrateNSteps(myStepper, *this, myY, snColl, 0.0, 0.0035, 10000);    
//   }
// };




  // bool success = true;
  // std::cout << std::boolalpha << success;

  // MPI_Init(&argc,&argv);
  // int rank; // My process ID
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Epetra_MpiComm Comm(MPI_COMM_WORLD);

  // int MyPID = Comm.MyPID();
  // int NumProc = Comm.NumProc();

  // int NumMyElements = 10000;
  // int NumMyElements1 = NumMyElements; // Needed for localmap
  // int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  // if (MyPID < 3) NumMyElements++;
  // int IndexBase = 0;
  // int ElementSize = 7;

  // Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase,Comm);
  // Epetra_Vector A(*LocalMap);

  // // epetramock::evector obj;
  // //obj.print();
 
  //   MPI_Finalize();
