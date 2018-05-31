#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>
// #include "ode_euler_stepper.hpp"
// #include "ode_rk4_stepper.hpp"
// #include "ode_integrate_n_steps.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

#include "apps_taygreen_epetra.hpp"


int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int rank = Comm.MyPID();
  //  MPI_Rank(MPI_COMM_WORLD, &rank);
  
  unsigned int Nx = 10;
  apps::taygreen app1(Comm, Nx);
  app1.run();
  //double err1 = app1.getError();
  
  // Nx = 400;
  // apps::taygreen app2(Comm, Nx);
  // app2.run();
  // double err2 = app2.getError();

  // if(rank==0)
  //   std::cout << std::log(err1/err2) << std::endl;
  
  MPI_Finalize() ;  
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
