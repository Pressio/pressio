
#include <gtest/gtest.h>
#include "Epetra_MpiComm.h"
#include "CORE_VECTOR"
#include "CORE_MATRIX"

// I need to finish this, for vectors it is easy to
// use maps to extract parts of it, but for crs matrix
// we need to be careful because we might need to use
// both column and row maps to extract parts of it.


TEST(core_trilinos_ma_export, Test1)
{
  //  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  auto rank = Comm.MyPID();
  int MyPID = Comm.MyPID();
  assert(Comm.NumProc() == 3);

  //CRS matrix
  int numGlobalEntries_ = 40;
  Epetra_Map Map(numGlobalEntries_, 0, Comm);
  Epetra_CrsMatrix A(Copy, Map, 1);
  auto myN1 = Map.NumMyElements();
  std::vector<int> mygel1(myN1);
  Map.MyGlobalElements(mygel1.data());
  std::array<double,1> vals({1});
  std::array<int,1> colind;
  for (auto & it : mygel1){
    colind[0] = it;
    vals[0] = it+100;
    A.InsertGlobalValues(it, 1, vals.data(), colind.data());
  }
  A.FillComplete();
  if (MyPID==0){
  std::cout << "=====================";
  std::cout << " AA ";
  std::cout << "=====================\n";}
  A.Print(std::cout);
  
  /////////////////////////////////////
  // target map
  int totN = 8;
  int myN = totN/Comm.NumProc();
  std::vector<int> mygel(myN);
  if(rank==0)
    mygel={0,4};
  if(rank==1)
    mygel={13,15};
  if(rank==2)
    mygel={21,25};
  if(rank==3)
    mygel={33,38};
  Epetra_Map destMap(-1, myN, mygel.data(), 0, Comm);
  //destMap.Print(std::cout);

  Map.Print(std::cout);
  destMap.Print(std::cout);
  
  Epetra_CrsMatrix B(Copy, destMap, 1);
  // // define importer: Epetra_Import(targetMap, sourceMap)
  Epetra_Import myImporter(destMap, Map);

  // // import global -> local
  B.Import(A, myImporter, Insert);
  B.FillComplete(Map, destMap);
  // B.OptimizeStorage();
  // std::cout << B.NumGlobalRows() << " " << B.NumGlobalCols() << std::endl;
  if (MyPID==0){
  std::cout << "=====================";
  std::cout << " BB ";
  std::cout << "=====================\n";}
  B.Print(std::cout);

}



TEST(core_trilinos_ma_export, Test2)
{
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  auto rank = Comm.MyPID();
  assert(Comm.NumProc() == 3);
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  int NumMyEquations = 5;
  int NumGlobalEquations = NumMyEquations*NumProc+3;
  if (MyPID < 3) NumMyEquations++;

  Epetra_Map SourceMap(NumGlobalEquations, NumMyEquations, 0, Comm);

  // Get update list and number of local equations from newly created Map
  int NumMyElements = SourceMap.NumMyElements();
  int * SourceMyGlobalElements = new int[NumMyElements];
  SourceMap.MyGlobalElements(SourceMyGlobalElements);
  // Construct a Target Map that will contain:
  //  some unchanged elements (relative to the soure map),
  //  some permuted elements
  //  some off-processor elements
  Epetra_Vector RandVec(SourceMap);
  RandVec.Random(); // This creates a vector of random numbers between negative one and one.

  int *TargetMyGlobalElements = new int[NumMyElements];
  int MinGID = SourceMap.MinMyGID();
  for (int i=0; i< NumMyEquations/2; i++) TargetMyGlobalElements[i] = i + MinGID; // Half will be the same...
  for (int i=NumMyEquations/2; i<NumMyEquations; i++) {
    int index = abs((int)(((double) (NumGlobalEquations-1) ) * RandVec[i]));
    TargetMyGlobalElements[i] = EPETRA_MIN(NumGlobalEquations-1,
					   EPETRA_MAX(0,index));
  }

  int NumSameIDs = 0;
  int NumPermutedIDs = 0;
  int NumRemoteIDs = 0;
  bool StillContiguous = true;
  for (int i=0; i < NumMyEquations; i++) {
    if (SourceMyGlobalElements[i]==TargetMyGlobalElements[i] && StillContiguous)
      NumSameIDs++;
    else if (SourceMap.MyGID(TargetMyGlobalElements[i])) {
      StillContiguous = false;
      NumPermutedIDs++;
    }
    else {
      StillContiguous = false;
      NumRemoteIDs++;
    }
  }
  Epetra_Map TargetMap(-1, NumMyElements, TargetMyGlobalElements, 0, Comm);

  if (MyPID==0){
  std::cout << "=====================";
  std::cout << " SOURCE MAP";
  std::cout << "=====================\n";}
  SourceMap.Print(std::cout);

  if (MyPID==0){
  std::cout << "=====================";
  std::cout << " TARGET MAP";
  std::cout << "=====================\n";}
  TargetMap.Print(std::cout);
  int i=0,j=0;;
  
  int NumVectors = 3;
  Epetra_MultiVector SourceMultiVector(SourceMap, NumVectors);
  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++)
      SourceMultiVector[j][i] = (double) SourceMyGlobalElements[i]*(j+1);

  // Create a target multivector that we will fill using an Import
  Epetra_MultiVector TargetMultiVector(TargetMap, NumVectors);
  Epetra_Import Importer(TargetMap, SourceMap);
  TargetMultiVector.Import(SourceMultiVector, Importer, Insert);

  if (MyPID==0){
    std::cout << "=====================";
    std::cout << " SOURCE VECTOR";
    std::cout << "=====================\n";
  }
  SourceMultiVector.Print(std::cout);

  if (MyPID==0){
  std::cout << "=====================";
  std::cout << " TARGET VECTOR";
  std::cout << "=====================\n";}
  TargetMultiVector.Print(std::cout);
  
  // // Now use Importer to do an export
  // Epetra_Vector TargetVector(SourceMap);
  // Epetra_Vector ExpectedTarget(SourceMap);
  // Epetra_Vector SourceVector(TargetMap);

  // NumSameIDs = Importer.NumSameIDs();
  // int NumPermuteIDs = Importer.NumPermuteIDs();
  // int NumExportIDs = Importer.NumExportIDs();
  // int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  // int *ExportLIDs = Importer.ExportLIDs();
  // int *ExportPIDs = Importer.ExportPIDs();
  // for (i=0; i < NumSameIDs; i++) ExpectedTarget[i] = (double) (MyPID+1);
  // for (i=0; i < NumPermuteIDs; i++) ExpectedTarget[PermuteFromLIDs[i]] =
  // 				      (double) (MyPID+1);
  // for (i=0; i < NumExportIDs; i++) ExpectedTarget[ExportLIDs[i]] +=
  // 				     (double) (ExportPIDs[i]+1);

  // for (i=0; i < NumMyElements; i++) SourceVector[i] =  (double) (MyPID+1);

  // TargetVector.Export(SourceVector, Importer, Insert);

  // if (MyPID==0){
  //   std::cout << "=====================";
  //   std::cout << " SOURCE VECTOR";
  //   std::cout << "=====================\n";
  // }
  // SourceVector.Print(std::cout);

  // if (MyPID==0){
  // std::cout << "=====================";
  // std::cout << " TARGET VECTOR";
  // std::cout << "=====================\n";}
  // TargetVector.Print(std::cout);

  
  
}
