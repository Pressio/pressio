
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>

class burg1d
{
private:
  using vecD = std::vector<double>;
  const double xL_ = 0.0;
  const double xR_ = 100.0;
  const vecD mu_ {5.0, 0.02, 0.02};
  const unsigned int nel_ = 10;
  const double t0 = 0.0;
  const double tfinal_ = 35.0;
  double dx_;
  vecD xGrid_;
  vecD U_;

public:
  using state_type = vecD;
  
  burg1d(){
    dx_ = (xR_ - xL_)/static_cast<double>(nel_);
    xGrid_.resize(nel_,0.0);
    for (int i=0; i<nel_; ++i){
      xGrid_[i] = dx_*i + dx_*0.5;
    };
  }

  void init(){
    U_.resize(nel_, 1.0);
  };

  void operator() ( const state_type & u,
		    state_type & R,
		    const double /* t */ )
  {
    R[0] = (0.5 * ( mu_[0]*mu_[0] - u[0]*u[0] ) )/dx_;
    for (int i=1; i<nel_; ++i){
      R[i] = ( 0.5*(u[i-1]*u[i-1] - u[i]*u[i]) )/dx_;
    }
    for (int i=0; i<nel_; ++i){
      R[i] += mu_[1]*exp(mu_[2]*xGrid_[i]);
    }    
  }
  void integrate()
  {
    int nsteps = 2500;
    vecD rhs(nel_,0.0);
    double dt = tfinal_/static_cast<double>(nsteps);

    std::ofstream file;
    file.open( "out.txt" );
    for (int step=0; step<nsteps; ++step)
    {
      (*this)(U_, rhs, step*dt);
      for (int i=0; i<nel_; ++i){
      	U_[i] += dt*(rhs[i]);
	if (step % 50 == 0 || step==0)
	 file << std::fixed << std::setprecision(10) << U_[i] << " ";
      }      
      if (step % 50 == 0)
	file << std::endl;
    }
    file.close();
  }

};


int main(int argc, char *argv[])
{
  burg1d app;
  app.init();
  app.integrate();
  
  return 0;
}


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
