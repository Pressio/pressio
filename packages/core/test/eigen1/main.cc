
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "core_ConfigDefs.hpp"
#include "meta.hpp"
#include "forwardDeclarations.hpp"
#include "staticAssert.hpp"

#include "vector/core_vectorTraits.hpp"
#include "vector/core_vector_epetra.hpp"
#include "vector/core_vector_serarb.hpp"
#include "vector/core_vector_std.hpp"

#include "matrix/core_matrixTraits.hpp"
#include "matrix/core_matrix_eigen.hpp"
// #include "matrix/core_matrix_std.hpp"

#include <Eigen/Dense>

template<typename T>
using vecT = std::vector<T>;

//using Eigen::MatrixXd;
int main()
{
  Eigen::MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;

  using mat_t = Eigen::Matrix<double,4,5>;
  mat_t AA;
  core::matrix<mat_t> BB;

  // using wrap_t = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  // using mat_t = core::matrix<wrap_t,double>;
  // mat_t A;

  // using wrap_t2 = vecT<vecT<double>>;
  // using mat_t2 = core::matrix<wrap_t2,double>;
  // mat_t2 A2;
  
  // using wrap_t3 = Eigen::Matrix<double,3,3>;
  // core::matrix<wrap_t3,double,3,3> DD;
  // //mat_t3 A3;

  
  return 0;
};




// // Epetra headers
// #include "Epetra_Time.h"
// #include "Epetra_BlockMap.h"
// #include "Epetra_Map.h"
// #include "Epetra_LocalMap.h"
// #include "Epetra_Vector.h"
// #include "Epetra_Version.h"
// #include "Epetra_MpiComm.h"
// #include <mpi.h>
// #include "Teuchos_VerboseObject.hpp"


// class myve{
// public:
//  double data_[5];
// public:
//   myve(){
//     data_[0] = 2.2;
//     data_[1] = 2.2;
//   };
//   ~myve(){};

//   double dot(const myve &) const{
//     return 1;
//   };
//   size_t size() const {
//     return 5;
//   };
//   double & operator [] (int i){
//     return data_[i];
//   };
//   double const & operator [] (int i) const{
//     return data_[i];
//   };
//   void resize(size_t val) {
//   }
// };


// int main(int argc, char *argv[]){
//   {
//     using state_t = std::vector<double>;
//     using vec_t = core::vector<state_t,double,int>;
//     vec_t gigi;
//     gigi.resize(10);
//     gigi[0] = 22.4;
//     state_t const * myview = gigi.view();
//     std::cout << (*myview)[0] << " " << (*myview)[1] << std::endl;
//     std::cout << core::details::traits<vec_t>::isSTDVector << std::endl;
//   }

//   {
//     using state_t = myve;
//     using vec_t = core::vector<state_t,double,int>;
//     vec_t gigi2;
//     size_t ss = gigi2.size();
//     std::cout << "size l = " << ss << std::endl;
//     auto * vv = gigi2.view();
//     std::cout << "view = " << (*vv)[0] << std::endl;
//   }

  
//   {  
//     MPI_Init(&argc,&argv);
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     Epetra_MpiComm Comm(MPI_COMM_WORLD);
//     int MyPID = Comm.MyPID();
//     int NumProc = Comm.NumProc();
//{
//     using lo_t = int;
//     using go_t = int;
    
//     const go_t numGlobalEntries = Comm.NumProc () * 10;
//     Epetra_Map contigMap (numGlobalEntries, 0, Comm);
  
//     Epetra_Vector x (contigMap);
//     x.PutScalar(21.0);


//     using vec_t = core::vector<Epetra_Vector,double,lo_t,go_t,
// 			       Epetra_Map,Epetra_MpiComm>;
//     vec_t gigi2(x);
//     std::cout << "size l = " << gigi2.localSize() << std::endl;
//     std::cout << "size g = " << gigi2.globalSize() << std::endl;
//     // // auto * vv = gigi2.view();
//     // // std::cout << "view = " << (*vv)[0] << std::endl;

//     // // Define some constants for use below.
//     // const double alpha = 3.14159;
//     // const double beta = 2.71828;
//     // const double gamma = -10.0;
//     // // x = beta*x + alpha*z
//     // x.Update (alpha, z, beta);
//     // y.PutScalar (42.0); // Set all entries of y to 42.0
//     // // y = gamma*y + alpha*x + beta*z
//     // y.Update (alpha, x, beta, z, gamma);
//     // // Compute the 2-norm of y.
//     // double theNorm = 0.0;
//     // y.Norm2 (&theNorm);

    
//     MPI_Finalize();
//   }
  
//   return 0;
// }
