
#include "laplace1dEpetra.hpp"

std::shared_ptr<Epetra_CrsMatrix> Laplace1dEpetra::calculateLinearSystem(){
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


std::shared_ptr<Epetra_Vector> Laplace1dEpetra::calculateForcingTerm(){
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

std::shared_ptr<Epetra_Vector> Laplace1dEpetra::createStates(){
  //---------------------------------------------------------------------------
  // Empty state Epetra_Vector for solver
  //---------------------------------------------------------------------------
  state_u = std::make_shared<nativeVec>(*contigMap);
  return state_u;
}

void Laplace1dEpetra::solveForStates(rcp<nativeMatrix> A, rcp<nativeVec> u, 
				     rcp<nativeVec> f){
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


void Laplace1dEpetra::printStates(rcp<nativeVec> u){
  //---------------------------------------------------------------------------
  // Print x and states                                                     
  //---------------------------------------------------------------------------
  for (int i = 0; i<nodesPerProc; i++){                            
    x_i = domain_[0]+domain_[2]*MyGlobalNodes[i];                        
    std::cout<<"x("<<MyGlobalNodes[i]<<"):\t" << x_i <<"\t"                    
             <<"u("<<MyGlobalNodes[i]<<"):\t" << (*u)[i] <<std::endl;
  }
}

std::shared_ptr<Epetra_Vector> Laplace1dEpetra::calcManufacturedForcing(){
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

void Laplace1dEpetra::compare2manufacturedStates(rcp<nativeVec> uapprox){
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

double Laplace1dEpetra::verifyImplementation(rcp<nativeMatrix> A){
  double norm = 0;
  std::shared_ptr<Epetra_Vector> fManu;
  std::shared_ptr<Epetra_Vector> uManuApprox;
  fManu = Laplace1dEpetra::calcManufacturedForcing(); 
  uManuApprox = Laplace1dEpetra::createStates();
  solveForStates(A, uManuApprox, fManu);
  Laplace1dEpetra::compare2manufacturedStates(uManuApprox);
  //---------------------------------------------------------------------------
  // Print the error between calculated states and manufactured solutions
  //---------------------------------------------------------------------------
  (void) (*uManuApprox).Norm2(&norm);
  return norm;
}
