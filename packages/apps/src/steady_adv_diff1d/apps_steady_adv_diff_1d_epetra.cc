
#include "apps_steady_adv_diff_1d_epetra.hpp"

namespace rompp{ namespace apps{

std::shared_ptr<Epetra_CrsMatrix> SteadyAdvDiff1dEpetra::calculateLinearSystem(){

  /* Couple comments:
   * (a) this method is doing too many things. Break it up more.
   * For instance, the creation of the map is something that not necessarily
   * relates to the linear system, so it should go in its own method.
   * To be fully precise, what you do below should be broken into few functions:
   *	(1) createMap
   *	(2) createMatrix
   *	(3) fillMatrix
   * the more you break up things, the better. For maintaining it, for clarity, etc.

   * (b) please avoid using raw pointers to manager memory. In c++11 you should
   * use shared_ptr or unique_ptr. In this case it is not too bad and it would not
   * create issues, but please get into the habit of NOT using new/delete.

   * (c) I don't know why calling "FillComplete" is a formality :)
   * Actually, calling fillComplete is very important.

   * (d) please rememeber to always initialize to zero variables when you declare them.
   * in C++11 this is easy, you can simple use {}.
   * for instance, say you have a int variable called a, you can do:
   * int a {};   or int a = {};   //this initializes it to 0
   * sample applies to objects of classes
   */

  //---------------------------------------------------------------------------
  //Create map
  //---------------------------------------------------------------------------
  numGlobalNodes = (domain_[1]-domain_[0])/domain_[2] + 1;
  //Make the contiguous Epetra Map with base of 0
  contigMap = std::make_shared<Epetra_Map>(numGlobalNodes, 0, *comm_);
  //Double check that the map is actually contiguous
  if (!contigMap->LinearMap())
    throw std::logic_error("The supposedly contiguous map isn't so.\n");
  //---------------------------------------------------------------------------
  //Find local nodes information
  //---------------------------------------------------------------------------
  nodesPerProc= contigMap->NumMyElements();
  MyGlobalNodes = contigMap->MyGlobalElements();
  NumNz = new int[nodesPerProc];      //empty to track how many non zero values
  for (int i=0; i<nodesPerProc; i++)
    if (MyGlobalNodes[i] ==0 || MyGlobalNodes[i]== numGlobalNodes-1)
      NumNz[i] = 2;
    else
      NumNz[i] = 3;
  //---------------------------------------------------------------------------
  //Create the sparse matrix A
  //---------------------------------------------------------------------------
  //Define individual values
  one = 1.0;                          //Boundary conditions
  two = (-mu_[0]*2.0)/(domain_[2]*domain_[2]); //Main diagonal
  Values = new double[2];             //Off diagonal values
  Indices = new int[2];
  Values[0] = (mu_[0]*1.0)/(domain_[2]*domain_[2])-(mu_[1]*1)/(2*domain_[2]);
  Values[1] = (mu_[0]*1.0)/(domain_[2]*domain_[2])+(mu_[1]*1)/(2*domain_[2]);

  //Create A and fill it in
  A = std::make_shared<nativeMatrix>(Copy, *contigMap, NumNz);

  //Iterate through number each local node/rows
  for (int i=0; i<nodesPerProc; i++){  //Set up the off diagonals
    if (MyGlobalNodes[i]> 0 && MyGlobalNodes[i]<numGlobalNodes-1){
      Indices[0] = MyGlobalNodes[i] - 1;  //Global Indicies
      Indices[1] = MyGlobalNodes[i] + 1;  //Global Indicies
      NumEntries = 2;
      //Fill in the off diagonals
      A->InsertGlobalValues(MyGlobalNodes[i], NumEntries, Values, Indices);
      //Indicies refer to the column index number here   iagonal 2
      A->InsertGlobalValues(MyGlobalNodes[i], 1, &two, MyGlobalNodes+i);
    }else{
      //Fill in the indicies that refer to the boundary conditions
      A->InsertGlobalValues(MyGlobalNodes[i], 1, &one, MyGlobalNodes+i);
    }
  }
  //Formalities
  A->FillComplete();
  //---------------------------------------------------------------------------
  //Free Memory
  //---------------------------------------------------------------------------
  delete(NumNz);
  delete(Values);
  delete(Indices);

  return A;
}


std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::calculateForcingTerm(){

  /* for instance, here, you call this method calculateForcingTerm
   * but you also create f inside. This is misleadding. You should only do
   * in here what you say you are doing. If you say calculate, then here you
   * should only calcualate and fill, and put the creation of f somewehre else.
   * Because f is a data member so it is visible everywhere inside the class.
   */

  //---------------------------------------------------------------------------
  //Create Epetra Vector for the forcing term
  //---------------------------------------------------------------------------
  f = std::make_shared<nativeVec>(*contigMap);
  //Iterate through all the nodes for each processor
  for (int i = 0; i<nodesPerProc; i++){
    if (MyGlobalNodes[i]> 0 && MyGlobalNodes[i]<numGlobalNodes-1){
      x_i = domain_[0]+domain_[2]*MyGlobalNodes[i]; //calculate x_i;
      (*f)[i] = exp(mu_[2]*x_i)*(3*x_i+x_i*x_i);    //has third parameter
    }else if (MyGlobalNodes[i]==0){
      (*f)[i] = bc1D_[0];                           //boundary condiitons Left
    }else if (MyGlobalNodes[i]==numGlobalNodes-1){
      (*f)[i] = bc1D_[1];                           //boundary conditions Right
    }
  }
  return f;
}

std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::createStates(){
  //---------------------------------------------------------------------------
  // Empty state Epetra_Vector for solver
  //---------------------------------------------------------------------------
  state_u = std::make_shared<nativeVec>(*contigMap);
  return state_u;
}

void SteadyAdvDiff1dEpetra::solveForStates(rcp<nativeMatrix> A, rcp<nativeVec> u,
				     rcp<nativeVec> f){
  /* options and params are not used, so remove them if not used */
  int options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];

  //---------------------------------------------------------------------------
  // Set, Create, and Solve Linear System
  //---------------------------------------------------------------------------
  Epetra_LinearProblem Problem(A.get(), u.get(), f.get());
  //Create AztecOO object
  AztecOO Solver(Problem);
  //AZ_defaults(options, params);
  //  Solver.SetAztecOption(AZ_precond, AZ_lu);
  //Solver.SetAllAztecParams(params);
  Solver.Iterate(100, 1E-6);
  Solver.NumIters();
  Solver.TrueResidual();
}


void SteadyAdvDiff1dEpetra::printStates(rcp<nativeVec> u){
  //---------------------------------------------------------------------------
  // Print x and states
  //---------------------------------------------------------------------------
  for (int i = 0; i<nodesPerProc; i++){
    x_i = domain_[0]+domain_[2]*MyGlobalNodes[i];
    std::cout<<"x("<<MyGlobalNodes[i]<<"):\t" << x_i <<"\t"
             <<"u("<<MyGlobalNodes[i]<<"):\t" << (*u)[i] <<std::endl;
  }
}

std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::calcManufacturedForcing(){
  //---------------------------------------------------------------------------
  //Create Epetra Vector for the forcing term
  //---------------------------------------------------------------------------
  f = std::make_shared<nativeVec>(*contigMap);
  double dLen = domain_[1]-domain_[0];
  //Iterate through all the nodes for each processor
  for (int i = 0; i<nodesPerProc; i++){
    if (MyGlobalNodes[i]> 0 && MyGlobalNodes[i]<numGlobalNodes-1){
      x_i = domain_[0]+domain_[2]*MyGlobalNodes[i]; //calculate x_i;
      (*f)[i] =
	mu_[0]*4*pow(M_PI,2)/pow(dLen,2)*cos(2*M_PI*(x_i-domain_[0])/dLen)+
	mu_[1]*2.0*M_PI/dLen*sin(2*M_PI*(x_i-domain_[0])/dLen);
    }else if (MyGlobalNodes[i]==0){
      (*f)[i] = 0;                           //boundary condiitons Left
    }else if (MyGlobalNodes[i]==numGlobalNodes-1){
      (*f)[i] = 0;                           //boundary conditions Right
    }
  }
  return f;
}

void SteadyAdvDiff1dEpetra::compare2manufacturedStates(rcp<nativeVec> uapprox){
  //---------------------------------------------------------------------------
  //Compare the manufactured solution to the discretized manufactured solution
  //---------------------------------------------------------------------------
  double uManu_i;
  for(int i=0; i<nodesPerProc; i++){
    x_i = domain_[0]+domain_[2]*MyGlobalNodes[i]; //calculate x_i;
    uManu_i = 1-cos(2*M_PI*(x_i-domain_[0])/(domain_[1]-domain_[0]));
    if (x_i> domain_[0] && x_i < domain_[1]){
      (*uapprox)[i] = fabs((*uapprox)[i]-uManu_i)/uManu_i;
    }else{
      (*uapprox)[i] = 0;
    }
  }
}

double SteadyAdvDiff1dEpetra::verifyImplementation(rcp<nativeMatrix> A){
  double norm = 0;
  std::shared_ptr<Epetra_Vector> fManu;
  std::shared_ptr<Epetra_Vector> uManuApprox;

  /* you don't need to do SteadyAdvDiff1dEpetra::
  * because you are inside the class, so eveyrthing is visibile.
  */
  fManu = SteadyAdvDiff1dEpetra::calcManufacturedForcing();
  uManuApprox = SteadyAdvDiff1dEpetra::createStates();
  solveForStates(A, uManuApprox, fManu);
  SteadyAdvDiff1dEpetra::compare2manufacturedStates(uManuApprox);

  //---------------------------------------------------------------------------
  // Print the error between calculated states and manufactured solutions
  //---------------------------------------------------------------------------
  /* why this? (void)? you can simply do:
     uManuApprox->Norm2(&norm);
  */
  (void) (*uManuApprox).Norm2(&norm);
  return norm;
}

}} //namespace rompp::apps
