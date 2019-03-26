
#include "apps_steady_adv_diff_1d_epetra.hpp"

namespace rompp{ 
  namespace apps{
    
    void SteadyAdvDiff1dEpetra::createMap(){
      //-----------------------------------------------------------------------
      //Create map to be used for the matrices and vectors 
      //-----------------------------------------------------------------------
      //Calculate the total number of global nodes
      numGlobalNodes_ = (domain_[1]-domain_[0])/domain_[2] + 1;
      //Make the contiguous Epetra Map with base of 0
      contigMap_ = std::make_shared<Epetra_Map>(numGlobalNodes_, 0, *comm_);
      //Double check that the map is actually contiguous
      if (!contigMap_->LinearMap())
	throw std::logic_error("The supposedly contiguous map isn't so.\n");
    }
    
    void SteadyAdvDiff1dEpetra::setup(){
      //-----------------------------------------------------------------------
      // Create Map and define indices of nonzero regions for linear system
      //-----------------------------------------------------------------------
      createMap();      
      nodesPerProc_= contigMap_->NumMyElements();     //Number of local nodes
      MyGlobalNodes_ = contigMap_->MyGlobalElements();//Retrives the global node
      std::unique_ptr<int[]> NumNz(new int[numGlobalNodes_]);
      for (int i=0; i<nodesPerProc_; i++)
      	if (MyGlobalNodes_[i] ==0 || MyGlobalNodes_[i]== numGlobalNodes_-1)
      	  NumNz[i] = 1;                             //Boundary conditions
      	else
      	  NumNz[i] = 3;                             //Inner nodes
      //-----------------------------------------------------------------------
      // Create the matrix and vectors to later fill
      //-----------------------------------------------------------------------
      //Create forcing term
      f_ = std::make_shared<nativeVec>(*contigMap_);
      //Create A
      A_ = std::make_shared<nativeMatrix>(Copy, *contigMap_, NumNz.get());
      //Create states and leave empty
      u_ = std::make_shared<nativeVec>(*contigMap_);
    }
    
    std::shared_ptr<Epetra_CrsMatrix> 
    SteadyAdvDiff1dEpetra::calculateLinearSystem(){ 
      //-----------------------------------------------------------------------
      //Calculate the components of A
      //-----------------------------------------------------------------------
      //Define individual values
      int NumEntries = 2;                                 //Off diagonal values
      double one = 1.0;                                   //Boundary conditions
      double two = (-mu_[0]*2.0)/(domain_[2]*domain_[2]); //Main diagonal
      //Define the diagonal values for advection and diffusion
      std::unique_ptr<double[]> Values(new double[2]);
      std::unique_ptr<int[]> Indices(new int[2]);
      Values[0] = (mu_[0]*1.0)/(domain_[2]*domain_[2])
	-(mu_[1]*1)/(2*domain_[2]);
      Values[1] = (mu_[0]*1.0)/(domain_[2]*domain_[2])
	+(mu_[1]*1)/(2*domain_[2]);
      //-----------------------------------------------------------------------
      //Fill in the values of Matrix A
      //-----------------------------------------------------------------------
      for (int i=0; i<nodesPerProc_; i++){  //Set up the off diagonals
	if (MyGlobalNodes_[i]> 0 && MyGlobalNodes_[i]<numGlobalNodes_-1){
	  Indices[0] = MyGlobalNodes_[i] - 1;  //Global Indicies
	  Indices[1] = MyGlobalNodes_[i] + 1;  //Global Indicies
	  //Fill in the off diagonals
	  A_->InsertGlobalValues(MyGlobalNodes_[i], NumEntries, Values.get(), 
				Indices.get());
	  //Indicies refer to the column index number here   iagonal 2
	  A_->InsertGlobalValues(MyGlobalNodes_[i], 1, &two, MyGlobalNodes_+i);
	}else{
	  //Fill in the indicies that refer to the boundary conditions
	  A_->InsertGlobalValues(MyGlobalNodes_[i], 1, &one, MyGlobalNodes_+i);
	}
      }
      //Finalize A Matrix Structure
      A_->FillComplete();
      
      return A_;
    }

    std::shared_ptr<Epetra_Vector> 
    SteadyAdvDiff1dEpetra::calculateForcingTerm(){  
      //-----------------------------------------------------------------------
      //Create Epetra Vector for the forcing term
      //-----------------------------------------------------------------------
      //Iterate through all the nodes for each processor
      for (int i = 0; i<nodesPerProc_; i++){
	if (MyGlobalNodes_[i]> 0 && MyGlobalNodes_[i]<numGlobalNodes_-1){
	  x_i_ = domain_[0]+domain_[2]*MyGlobalNodes_[i]; //calculate x_i;
	  (*f_)[i] = exp(mu_[2]*x_i_)*(3*x_i_+x_i_*x_i_); //has third parameter
	}else if (MyGlobalNodes_[i]==0){
	  (*f_)[i] = bc1D_[0];                      //boundary condiitons Left
	}else if (MyGlobalNodes_[i]==numGlobalNodes_-1){
	  (*f_)[i] = bc1D_[1];                      //boundary conditions Right
	}
      }
      
      return f_;
    }
    
    std::shared_ptr<Epetra_Vector> SteadyAdvDiff1dEpetra::getStates(){
      //-----------------------------------------------------------------------
      // Empty state Epetra_Vector for solver
      //-----------------------------------------------------------------------
      return u_;
    }
    
    void SteadyAdvDiff1dEpetra::calculateStates(rcp<nativeMatrix> A, 
					       rcp<nativeVec> u, 
					       rcp<nativeVec> f){
      //-----------------------------------------------------------------------
      // Set, Create, and Solve Linear System
      //-----------------------------------------------------------------------
      Epetra_LinearProblem Problem(A.get(), u.get(), f.get());
      //Create AztecOO object
      AztecOO Solver(Problem);
      Solver.Iterate(100, 1E-6);
      Solver.NumIters(); 
      Solver.TrueResidual();
    }
    
    void SteadyAdvDiff1dEpetra::printStates(rcp<nativeVec> u){
      //-----------------------------------------------------------------------
      // Print x and states
      //-----------------------------------------------------------------------
      for (int i = 0; i<nodesPerProc_; i++){
	x_i_ = domain_[0]+domain_[2]*MyGlobalNodes_[i];
	std::cout<<"x("<<MyGlobalNodes_[i]<<"):\t" << x_i_ <<"\t"
		 <<"u("<<MyGlobalNodes_[i]<<"):\t" << (*u)[i] <<std::endl;
      }
    }
    
    std::shared_ptr<Epetra_Vector> 
    SteadyAdvDiff1dEpetra::calculateManufacturedForcing(){
      //-----------------------------------------------------------------------
      //Create Epetra Vector for the forcing term
      //-----------------------------------------------------------------------
      f_ = std::make_shared<nativeVec>(*contigMap_);
      double dLen = domain_[1]-domain_[0];
      //Iterate through all the nodes for each processor
      for (int i = 0; i<nodesPerProc_; i++){
	if (MyGlobalNodes_[i]> 0 && MyGlobalNodes_[i]<numGlobalNodes_-1){
	  x_i_ = domain_[0]+domain_[2]*MyGlobalNodes_[i]; //calculate x_i;
	  (*f_)[i] = mu_[0]*4*pow(M_PI,2)/pow(dLen,2)*
	    cos(2*M_PI*(x_i_-domain_[0])/dLen)+
	    mu_[1]*2.0*M_PI/dLen*sin(2*M_PI*(x_i_-domain_[0])/dLen);
	}else if (MyGlobalNodes_[i]==0){
	  (*f_)[i] = 0;                           //boundary condiitons Left
	}else if (MyGlobalNodes_[i]==numGlobalNodes_-1){
	  (*f_)[i] = 0;                           //boundary conditions Right
	}
      }
      return f_;
    }
    
    void 
    SteadyAdvDiff1dEpetra::compare2manufacturedStates(rcp<nativeVec> uapprox){
      //-----------------------------------------------------------------------
      //Compare the manufactured state to the discretized manufactured solution
      //-----------------------------------------------------------------------
      double uManu_i;
      for(int i=0; i<nodesPerProc_; i++){
	x_i_ = domain_[0]+domain_[2]*MyGlobalNodes_[i]; //calculate x_i;
	uManu_i = 1-cos(2*M_PI*(x_i_-domain_[0])/(domain_[1]-domain_[0]));
	if (x_i_> domain_[0] && x_i_ < domain_[1]){
	  (*uapprox)[i] = fabs((*uapprox)[i]-uManu_i)/uManu_i;
	}else{
	  (*uapprox)[i] = 0;
	}
      }
    }
    
    double SteadyAdvDiff1dEpetra::verifyImplementation(rcp<nativeMatrix> A){
      //-----------------------------------------------------------------------
      // Verify that approximate states match with the manufactured solution
      //-----------------------------------------------------------------------
      double norm = 0;
      std::shared_ptr<Epetra_Vector> fManu;
      std::shared_ptr<Epetra_Vector> uManuApprox;
      
      fManu = calculateManufacturedForcing();  //Calculate manu forcing term
      uManuApprox = getStates();               //Get empty states
      calculateStates(A, uManuApprox, fManu);  //Solve for the manu states
      compare2manufacturedStates(uManuApprox); //Compare approx states to manu
      
      //-----------------------------------------------------------------------
      // Compare error between calculated states and manufactured solutions
      //-----------------------------------------------------------------------
      (*uManuApprox).Norm2(&norm);             // L2 norm of the error

      return norm;
    }    
  }
} //namespace rompp::apps
