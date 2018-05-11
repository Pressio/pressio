
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
// #include "BuildTestProblems.h"
// #include "ExecuteTestProblems.h"
// #include "../epetra_test_err.h"
#include "Epetra_Version.h"
// #ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#  include <mpi.h>
// #else
// #  include "Epetra_SerialComm.h"
// #endif

#  include "Teuchos_VerboseObject.hpp"


int main(int argc, char *argv[])
{
  bool success = true;
  std::cout << std::boolalpha << success;

  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();



  int NumMyElements = 10000;
  int NumMyElements1 = NumMyElements; // Needed for localmap
  int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  int ElementSize = 7;

  Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase,Comm);
  Epetra_Vector A(*LocalMap);

  // epetramock::evector obj;
  //obj.print();
 
   MPI_Finalize();

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
