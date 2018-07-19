
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "core_ConfigDefs.hpp"
#include "core_meta.hpp"
#include "core_forward_declarations.hpp"
#include "core_static_assert.hpp"

#include "vector/core_vector_traits.hpp"
#include "vector/core_vector_epetra.hpp"
#include "vector/core_vector_serial_arbitrary.hpp"
#include "vector/core_vector_std.hpp"
#include "vector/core_vector_eigen.hpp"

// Epetra headers
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_MpiComm.h"
#include <mpi.h>
#include "Teuchos_VerboseObject.hpp"

#include "Eigen/Dense"


class myve{
public:
 double data_[5];
public:
  using scalar_type = double;
  using ordinal_type = int;
  myve(){
    data_[0] = 2.2;
    data_[1] = 2.2;
  };
  ~myve(){};

  double dot(const myve &) const{
    return 1;
  };
  size_t size() const {
    return 5;
  };
  double & operator [] (int i){
    return data_[i];
  };
  double const & operator [] (int i) const{
    return data_[i];
  };
  void resize(size_t val) {
  }
};


int main(int argc, char *argv[]){
  {
    using state_t = std::vector<double>;
    using vec_t = core::Vector<state_t>;
    vec_t gigi;
    gigi.resize(10);
    gigi[0] = 22.4;
    state_t const * myview = gigi.view();
    std::cout << (*myview)[0] << " " << (*myview)[1] << std::endl;
    std::cout << core::details::traits<vec_t>::isSTDVector << std::endl;
  }

  {
    using state_t = myve;
    using vec_t = core::Vector<state_t>;
    vec_t gigi2;
    size_t ss = gigi2.size();
    std::cout << "size l = " << ss << std::endl;
    auto * vv = gigi2.view();
    std::cout << "view = " << (*vv)[0] << std::endl;
  }

  
  {  
    MPI_Init(&argc,&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    const int numGlobalEntries = Comm.NumProc () * 10;
    Epetra_Map contigMap (numGlobalEntries, 0, Comm);
  
    Epetra_Vector x (contigMap);
    x.PutScalar(21.0);

    using vec_t = core::Vector<Epetra_Vector>;
    vec_t gigi2(x);
    std::cout << "size l = " << gigi2.localSize() << std::endl;
    std::cout << "size g = " << gigi2.globalSize() << std::endl;
    auto * vv = gigi2.view();
    std::cout << "view = " << (*vv)[0] << std::endl;

    // // Define some constants for use below.
    // const double alpha = 3.14159;
    // const double beta = 2.71828;
    // const double gamma = -10.0;
    // // x = beta*x + alpha*z
    // x.Update (alpha, z, beta);
    // y.PutScalar (42.0); // Set all entries of y to 42.0
    // // y = gamma*y + alpha*x + beta*z
    // y.Update (alpha, x, beta, z, gamma);
    // // Compute the 2-norm of y.
    // double theNorm = 0.0;
    // y.Norm2 (&theNorm);    
    MPI_Finalize();
   }

  {
    using mat_t = Eigen::Matrix<double,1,5>;
    mat_t M1;
    core::Vector<mat_t> V1;
  }

  
  return 0;
}





/*
// #define TEST_FIND_SUBSTR_IN_STR(SUBSTR, STR) \
//   { \
//     const bool foundSubStr = ((STR).find(SUBSTR) != std::string::npos); \
//     std::cout << "Found \"" SUBSTR "\" ? " << foundSubStr << "\n"; \
//     if (!foundSubStr) success=false; \
//   } \
//   (void)(success)


// int main() {

//   bool success = true;
//   std::cout << std::boolalpha;

//   SimpleCxx::HelloWorld helloWorld;
//   std::ostringstream oss;
//   helloWorld.printHelloWorld(oss);
//   std::cout << oss.str();

//   TEST_FIND_SUBSTR_IN_STR("Hello World", oss.str());

// #ifdef HAVE_SIMPLECXX___INT64
//   TEST_FIND_SUBSTR_IN_STR("We have __int64", oss.str());
// #endif

// #ifdef HAVE_SIMPLECXX_DEBUG
//   TEST_FIND_SUBSTR_IN_STR("Debug is enabled", oss.str());
// #else
//   TEST_FIND_SUBSTR_IN_STR("Release is enabled", oss.str());
// #endif

//   TEST_FIND_SUBSTR_IN_STR("Sqr(3) = 9", oss.str());

// #ifdef HAVE_SIMPLECXX_SIMPLETPL
//   TEST_FIND_SUBSTR_IN_STR("Cube(3) = 27", oss.str());
// #endif

//   if (success) {
//     std::cout << "End Result: TEST PASSED\n";
//   }
//   else {
//     std::cout << "End Result: TEST FAILED\n";
//   }
//   return (success ? EXIT_SUCCESS : EXIT_FAILURE);

// }
*/
