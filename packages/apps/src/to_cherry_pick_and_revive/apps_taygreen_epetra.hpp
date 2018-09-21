#ifndef APPS_TAYGREEN_EPETRA_HPP_
#define APPS_TAYGREEN_EPETRA_HPP_

#include "apps_ConfigDefs.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include "AztecOO.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>
#include <memory>
// #include "MueLu.hpp"
// #include "MueLu_FactoryManager.hpp"
// #include "MueLu_DirectSolver.hpp"
// #include "MueLu_Hierarchy.hpp"
// #include "MueLu_SaPFactory.hpp"
// #include "MueLu_GenericRFactory.hpp"
// #include "MueLu_RAPFactory.hpp"
// #include "MueLu_TentativePFactory.hpp"
// #include "MueLu_Ifpack2Smoother.hpp"
// #include "MueLu_SmootherFactory.hpp"
// #include <MueLu_TrilinosSmoother.hpp> 
// #include <MueLu_EpetraOperator.hpp>
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

namespace apps{

class taygreen
{
private:
  using ui_t = unsigned int;
  const double PI = 3.14159265358979323846;
  using vecui = std::vector<ui_t>;
  using vecd = std::vector<double>;
  using veci = std::vector<int>;

  template<typename T>
  using rcp = std::shared_ptr<T>;

public:
  taygreen(Epetra_MpiComm comm, ui_t N)
    : ecomm_(comm), gN_(N)
  {
    // myID_ = ecomm_.MyPID();
    // nProc_ = ecomm_.NumProc();
    // nPi_ = std::sqrt<ui_t>(nProc_);
    // nPj_ = nPi_;
    // kToij(myIDi_, myIDj_, nPj_, myID_);
    // lN_ = static_cast<ui_t>(gN_)/static_cast<ui_t>(nPi_);
    // lNsq_ = lN_*lN_;
    // dx_ = L_/static_cast<ui_t>(gN_-1);
    // // figure out what elemets I own
    // totN_ = gN_*gN_;
    // setup();
  }//end constructor

  // default constructor
  ~taygreen(){}


};

#endif 



  
  ui_t ijTok(ui_t i, ui_t j, ui_t N){
    return i*N + j;
  }
  void kToij(ui_t & i, ui_t & j, ui_t N, ui_t k){
    j = kToj(N,k);
    i = kToi(N,k);
  }
  ui_t kToi(ui_t N, ui_t k){
    return k/N;
  }
  ui_t kToj(ui_t N, ui_t k){
    return k % N;
  }
  
  void createMaps()
  {
    ui_t NumMyElements;
    veci MyGlobElements;
    ui_t myj_l = myIDj_*lN_;
    ui_t myj_r = myj_l + lN_;
    ui_t myi_l = myIDi_*lN_;
    ui_t myi_r = myi_l + lN_;

    NumMyElements = 0;
    for(ui_t i1=myi_l; i1<myi_r; ++i1){
     for(ui_t j1=myj_l; j1<myj_r; ++j1){
       ui_t gid = ijTok(i1,j1,gN_);
       MyGlobElements.emplace_back(gid);
       NumMyElements++;
     }
    }
    MyGlobElemPoissonMap_.resize(NumMyElements);
    noGhostMap_ = std::make_shared<Epetra_Map>(-1, (int) NumMyElements,
					       MyGlobElements.data(),
					       0, ecomm_ ) ;
    noGhostMap_->MyGlobalElements( MyGlobElemPoissonMap_.data() );

    //store left ghost buffer
    ui_t j1 = myIDj_==0 ? gN_-2 : myj_l-1;
    for(ui_t i1=myi_l; i1<myi_r; ++i1){
       ui_t gid = ijTok(i1,j1,gN_);
       MyGlobElements.emplace_back(gid);
       NumMyElements++;
    }
    //store right ghost buffer
    j1 = myIDj_==nPj_-1 ? 1 : myj_r;
    for(ui_t i1=myi_l; i1<myi_r; ++i1){
       ui_t gid = ijTok(i1,j1,gN_);
       MyGlobElements.emplace_back(gid);
       NumMyElements++;
    }    
    //store top ghost buffer
    ui_t i1 = myIDi_==0 ? gN_-2 : myi_l-1;
    for(ui_t j1=myj_l; j1<myj_r; ++j1){
       ui_t gid = ijTok(i1,j1,gN_);
       MyGlobElements.emplace_back(gid);
       NumMyElements++;
    }
    //store bottom ghost buffer
    i1 = myIDi_==nPi_-1 ? 1 : myi_r;
    for(ui_t j1=myj_l; j1<myj_r; ++j1){
       ui_t gid = ijTok(i1,j1,gN_);
       MyGlobElements.emplace_back(gid);
       NumMyElements++;
    }

    MyGlobElemFieldMap_.resize(NumMyElements);
    withGhostMap_ = std::make_shared<Epetra_Map>(-1, (int) NumMyElements,
					     MyGlobElements.data(),
					     0, ecomm_ ) ;
    withGhostMap_->MyGlobalElements( MyGlobElemFieldMap_.data() );

  }//end create map

  void fillGrid()
  {
    int NumMyElements = withGhostMap_->NumMyElements();
    int k = 0;
    for (auto const & it : MyGlobElemPoissonMap_){
      double x = dx_ * kToj(gN_, it);
      double y = L_ - dx_ * kToi(gN_, it);
      (*xx_)[k] = x;
      (*yy_)[k] = y;
      k++;
    }    
  }

  double vorticity(double x, double y, double t, double q = 1.0)
  {
    // double x1 = 3.0*PI/4.0;
    // double y1 = PI;
    // double x2 = 5.0*PI/4.0;
    // double y2 = PI;
    // double d1 = (x-x1)*(x-x1)+(y-y1)*(y-y1);
    // double d2 = (x-x2)*(x-x2)+(y-y2)*(y-y2);
    // return std::exp(-PI*d1) + std::exp(-PI*d2);
    // if (std::abs( x-PI )<0.5 && std::abs( y-PI )<0.5)
    //   return 1.0;
    // else
    //   return 0.0;
    return 2.0*std::sin(q*x)*std::sin(q*y)*std::exp(-2.*q*q*t/Re_);
  }
  double vel_U(double x, double y, double t, double q = 1.0){
    return std::sin(q*x)*std::cos(q*y)*std::exp(-2.*q*q*t/Re_);
  }
  double vel_V(double x, double y, double t, double q = 1.0){
    return -std::cos(q*x)*std::sin(q*y)*std::exp(-2.*q*q*t/Re_);
  }
  double stream(double x, double y, double t, double q = 1.0){
    return std::sin(q*x)*std::sin(q*y)*std::exp(-2.*q*q*t/Re_);
  }

  void printField(std::string which, rcp<Epetra_Vector> F, int step)
  {
    std::ofstream file;
    std::string fileName = which +
      "_p"+std::to_string(myID_)+
      "_s"+std::to_string(step)+".txt";      
    file.open( fileName );
    for(ui_t i=0; i < lN_; i++){
      for(ui_t j=0; j < lN_; j++){
	int kk = ijTok(i, j, lN_);
	file << std::fixed << std::setprecision(10) << (*F)[kk] << " ";
      }
      file << std::endl;
    }
    file.close();    
  };

  
public:
  taygreen(Epetra_MpiComm comm, ui_t N)
    : ecomm_(comm), gN_(N)
  {
    myID_ = ecomm_.MyPID();
    nProc_ = ecomm_.NumProc();
    nPi_ = std::sqrt<ui_t>(nProc_);
    nPj_ = nPi_;
    kToij(myIDi_, myIDj_, nPj_, myID_);
    lN_ = static_cast<ui_t>(gN_)/static_cast<ui_t>(nPi_);
    lNsq_ = lN_*lN_;
    dx_ = L_/static_cast<ui_t>(gN_-1);
    // figure out what elemets I own
    totN_ = gN_*gN_;
    setup();
  }//end constructor

  // default constructor
  ~taygreen(){}


  void fillLapla()
  {
    veci MyGlobalElements( lNsq_ );
    noGhostMap_->MyGlobalElements( MyGlobalElements.data() );
    //    noGhostMap_->Print(std::cout);
    veci Indices(4);
    vecd ValueDiag(1);
    veci IndicesDiag(1);
    double diagVal = -4.0;
    double westVal = 1.0;
    double eastVal = westVal, northVal=westVal, southVal=westVal;
    vecd Values{westVal, eastVal, northVal, southVal};

    int westID, eastID, northID, southID, diagID;
    for (int li=0; li<lN_; li++){
     for (int lj=0; lj<lN_; lj++)
     {
       int lkk = ijTok(li, lj, lN_); // local counter 
       int thisGID = MyGlobalElements[lkk]; // global ID of current point
       int gi = kToi(gN_, thisGID);
       int gj = kToj(gN_, thisGID);

       if (gi>1 && gi<gN_-1 && gj>1 && gj<gN_-1)
	 {
	   westID = (gj==0) ? ijTok(gi, gN_-2, gN_) : thisGID-1;
	   eastID = (gj==gN_-1) ? ijTok(gi, 1, gN_) : thisGID+1;
	   northID = (gi==0) ? ijTok(gN_-2, gj, gN_) : thisGID-gN_;
	   southID = (gi==gN_-1) ? ijTok(1, gj, gN_) : thisGID+gN_;
	   Indices[0] = westID; 
	   Indices[1] = eastID;
	   Indices[2] = northID;
	   Indices[3] = southID;
	   //off diagonal
	   lapMat_->InsertGlobalValues(thisGID, 4, Values.data(), Indices.data() );
	   // add diagonal
	   ValueDiag[0] = diagVal;
	   IndicesDiag[0] = thisGID;
	   lapMat_->InsertGlobalValues( thisGID, 1, ValueDiag.data(), IndicesDiag.data() );
	 }
       else
	 {
	   ValueDiag[0] = 1.0;
	   IndicesDiag[0] = thisGID;
	   lapMat_->InsertGlobalValues( thisGID, 1, ValueDiag.data(), IndicesDiag.data() );	   
	 }
       
     }//end j
    }//end i

    lapMat_->FillComplete();
    //    lapMat_->Print(std::cout);
  }//edn method 

  
  void setup()
  {
    // create maps
    createMaps();
    basicToGhostImporter_ = std::make_shared<Epetra_Import>(*withGhostMap_, *noGhostMap_);
    ghostToBasicImporter_ = std::make_shared<Epetra_Import>(*noGhostMap_, *withGhostMap_);
    //noGhostMap_->Print(std::cout);

    // fill grid info 
    xx_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    yy_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    fillGrid();

    // initialize fields
    vortNoGhost_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    vortWithGhost_  = std::make_shared<Epetra_Vector>(*withGhostMap_);
    // vortWithGhost1_ = std::make_shared<Epetra_Vector>(*withGhostMap_);
    // vortWithGhost2_ = std::make_shared<Epetra_Vector>(*withGhostMap_);
    for (int i=0; i<lNsq_; i++){
      (*vortNoGhost_)[i] = vorticity( (*xx_)[i], (*yy_)[i], 0.0 );
      (*vortWithGhost_)[i] = (*vortNoGhost_)[i];
    }
    vortWithGhost_->Import(*vortNoGhost_, *basicToGhostImporter_, Insert);

    UUNoGhost_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    VVNoGhost_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    streamNoGhost_ = std::make_shared<Epetra_Vector>(*noGhostMap_);
    streamWithGhost_  = std::make_shared<Epetra_Vector>(*withGhostMap_);
   
    int numOfNonZerosEntries = 5;
    lapMat_ = std::make_shared<Epetra_CrsMatrix>(Copy, *noGhostMap_, numOfNonZerosEntries );
    fillLapla();
    
    // if(myID_==0){
    //   vortWithGhost_->Print(std::cout);
    // }
    // if(myID_==0){
    //   vortWithGhost_->Print(std::cout);
    //   // for (int i=0; i<NumMyElementsNoGhosts; i++){
    //   // 	std::cout << " i= " << i
    //   // 		  << " val= " << (*vortNoGhost_)[i]
    //   // 		  << std::endl;
    //   // }
    // }   
  }

  void calcRHS(vecd & G, rcp<Epetra_Vector> ome, double t)
  {
    double north, south, west, east;
    for (int i=0; i<lN_; i++){
     for (int j=0; j<lN_; j++)
     {
       int kk = ijTok(i, j, lN_);
       north = (i==0) ? (*ome)[lNsq_+2*lN_+j] : (*ome)[kk-lN_];
       south = (i==lN_-1) ? (*ome)[lNsq_+3*lN_+j] : (*ome)[kk+lN_];
       west = (j==0) ? (*ome)[lNsq_+i] : (*ome)[kk-1];
       east = (j==lN_-1) ? (*ome)[lNsq_+lN_+i] : (*ome)[kk+1];

       double UU = (*UUNoGhost_)[kk]; //vel_U( (*xx_)[kk], (*yy_)[kk], t );
       double VV = (*VVNoGhost_)[kk]; //vel_V( (*xx_)[kk], (*yy_)[kk], t );
       // double UU = vel_U( (*xx_)[kk], (*yy_)[kk], t );
       // double VV = vel_V( (*xx_)[kk], (*yy_)[kk], t );
       double convX = UU*( east-west )/(2.0*dx_); 
       double convY = VV*( north-south )/(2.0*dx_);
       double diffX = ReInv_*( east - 2*(*ome)[kk] + west )/(dx_*dx_); 
       double diffY = ReInv_*( north - 2*(*ome)[kk] + south )/(dx_*dx_); 
       G[kk] = -convX-convY + diffX + diffY;
     }
    }
  }


  void computeError(double time){
    double localres = 0.0;
    // each process computes local error 
    for (int i=0; i<lNsq_; i++){
      double err = (*vortWithGhost_)[i] - vorticity((*xx_)[i],(*yy_)[i],time);
      localres += err*err;
    }    
    l2error_ = 0.0;
    ecomm_.SumAll(&localres,&l2error_,1);
    l2error_ = std::sqrt( l2error_ * dx_ * dx_);
  }

  void solvePoisson()
  {
    for (int i=0; i<lNsq_; i++){
      (*vortNoGhost_)[i] *= -(dx_*dx_);
      (*streamNoGhost_)[i] = 0.0;//stream((*xx_)[i],(*yy_)[i],0.0);
    }      
    Epetra_LinearProblem Problem( &(*lapMat_), &(*streamNoGhost_), &(*vortNoGhost_) );
    AztecOO Solver(Problem);
    Solver.SetAztecOption( AZ_solver, AZ_GMRESR);
    Solver.SetAztecOption( AZ_conv, AZ_r0);
    //Solver.SetAztecOption( AZ_conv, AZ_noscaled);
    Solver.SetAztecOption( AZ_output, AZ_all );
    // Solver.SetAztecOption( AZ_precond, AZ_Jacobi );
    // Solver.SetAztecOption( AZ_poly_ord, 10 );

    // ================= MultiLevelOperator SECTION ========================
    int nLevels = 10;            // maximum number of levels
    int maxMgLevels = 6;         //
    ML_Set_PrintLevel(10);       // print level (0 silent, 10 verbose)
    ML* ml_handle;               // container of all ML' data
    ML_Create(&ml_handle, maxMgLevels);
    // convert to epetra matrix, put finest matrix into
    // position maxMgLevels - 1 of the hierarchy. NOTE: the matrix
    // is only wrapped (that is, a suitable getrow() function is used),
    // so data in the linear system matrix are NOT replicated.
    EpetraMatrix2MLMatrix(ml_handle, maxMgLevels-1, lapMat_.get() );
    // create an Aggregate object; this will contain information
    // about the aggregation process for each level
    ML_Aggregate *agg_object;
    ML_Aggregate_Create(&agg_object);
    // select coarsening scheme.
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(agg_object);
    // generate the hierarchy. We decided to use decreasing ordering;
    // one can also use ML_INCREASING (in this case, you need to replace
    // maxMgLevels-1 with 0 in call to EpetraMatrix2MLMatrix())
    nLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, maxMgLevels-1,
    						  ML_DECREASING, agg_object);
    // // define the ID of the coarsest level
    int coarsestLevel = maxMgLevels - nLevels;
    // // set up some smoothers. 
    int nits = 1;
    for (int level = maxMgLevels-1; level > coarsestLevel; level--){
    	ML_Gen_Smoother_Cheby(ml_handle, level, ML_BOTH, 30., 3);
     } 
    // simple coarse solver. You may want to use Amesos to access
    // to a large variety of direct solvers, serial and parallel
    ML_Gen_Smoother_GaussSeidel(ml_handle, coarsestLevel, ML_BOTH,
    				nits, ML_DEFAULT);
    // generate the solver
    ML_Gen_Solver(ml_handle, ML_MGV, maxMgLevels-1, coarsestLevel);
    // wrap ML_Operator into Epetra_Operator
    ML_Epetra::MultiLevelOperator  MLPrec(ml_handle, ecomm_, *noGhostMap_, *noGhostMap_);
    // ========== End of MultiLevelOperator SECTION ========================

    // set this operator as preconditioner for AztecOO
    Solver.SetPrecOperator(&MLPrec);

    // solve
    Solver.Iterate(5000, 1e-14);
    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_handle);
  
    //Solver.Iterate(5000, 1e-6);
    // vecd status( AZ_STATUS_SIZE );
    // Solver.GetAllAztecStatus( status.data() );
    // double AztecOOExitStatus = status.operator[](1);
    // std::cout << "AztecOOExitStatus = " << AztecOOExitStatus
    // 	      << " " << AZ_normal << std::endl;
  }


  void computeVelocity()
  {
    double north, south, west, east;
    for (int i=0; i<lN_; i++){
     for (int j=0; j<lN_; j++)
     {
       int kk = ijTok(i, j, lN_);
       north = (i==0) ? (*streamWithGhost_)[lNsq_+2*lN_+j] : (*streamWithGhost_)[kk-lN_];
       south = (i==lN_-1) ? (*streamWithGhost_)[lNsq_+3*lN_+j] : (*streamWithGhost_)[kk+lN_];
       west = (j==0) ? (*streamWithGhost_)[lNsq_+i] : (*streamWithGhost_)[kk-1];
       east = (j==lN_-1) ? (*streamWithGhost_)[lNsq_+lN_+i] : (*streamWithGhost_)[kk+1];
       (*UUNoGhost_)[kk] = ( north-south )/(2.0*dx_);
       (*VVNoGhost_)[kk] = -( east-west )/(2.0*dx_);	 
     }
    }      
  }
  
  void run()
  {
    vecd Gn(lNsq_, 0.0);

    printField("x",  xx_, 0);
    printField("y",  yy_, 0);
    printField("vort", vortWithGhost_, 0);

    solvePoisson();
    //    streamNoGhost_->Print(std::cout);
    for (int i=0; i<lNsq_; i++){
      std::cout << std::setprecision(14)
    		<< (*streamNoGhost_)[i]
    		<< " " << stream((*xx_)[i],(*yy_)[i],0.0)
    		<< std::endl;
    }
    
    //streamWithGhost_->Import(*streamNoGhost_, *basicToGhostImporter_, Insert);
    //    computeVelocity();

    // double time = 0.0;
    // for (ui_t step = 0; step < 15000; ++step)
    // {
    //   if (myID_==0){
    // 	std::cout << " ----------- " << std::endl;
    // 	std::cout << " Iteration = " << step << std::endl;
    //   }

    //   // RK 1
    //   calcRHS(Gn, vortWithGhost_, time);
    //   for (int i=0; i<lNsq_; i++){
    //   	(*vortWithGhost_)[i] = (*vortWithGhost_)[i] + dt_ * Gn[i];
    //   }

      // for (int i=0; i<lNsq_; i++){
      // 	(*vortNoGhost_)[i] = (*vortWithGhost_)[i];
      // }
    ////   vortNoGhost_->Import(*vortWithGhost_, *ghostToBasicImporter_, Insert);
    //   vortWithGhost_->Import(*vortNoGhost_, *basicToGhostImporter_, Insert);
    //   solvePoisson();
    //   streamWithGhost_->Import(*streamNoGhost_, *basicToGhostImporter_, Insert);
    //   computeVelocity();

    //   time+=dt_;

    //   if (step % 50 ==0){
    // 	printField("vort", vortWithGhost_, step);
    // 	//computeError(time);
    // 	// if (myID_==0)
    // 	//   std::cout << " step = " << step
    //   	// 	    << std::setprecision(14)
    // 	// 	    << " Error = " << l2error_
    // 	// 	    << std::endl;	
    //   }//end if
    // }//end for

  }//end run

  double getError(){
    return l2error_;
  }
  
private:
  // comm
  Epetra_MpiComm ecomm_;
  ui_t myID_;
  ui_t myIDi_, myIDj_;
  ui_t nProc_;
  ui_t nPi_, nPj_;
  // distribution  
  rcp<Epetra_Map> withGhostMap_;
  rcp<Epetra_Map> noGhostMap_;
  veci MyGlobElemFieldMap_;
  veci MyGlobElemPoissonMap_;
  rcp<Epetra_Import> basicToGhostImporter_;
  rcp<Epetra_Import> ghostToBasicImporter_;
  double l2error_;
  // domain info
  const double dt_ = 0.005;
  const double Re_ = 1000.;
  const double ReInv_ = 1.0/Re_;
  const double L_ = 2.*PI;
  const double x0 = 0.0;
  const double y0 = x0;
  // grid info 
  ui_t gN_;
  ui_t totN_;
  ui_t lN_;
  ui_t lNsq_;
  double dx_;
  // vorticity fields
  rcp<Epetra_Vector> vortNoGhost_;
  rcp<Epetra_Vector> vortWithGhost_;
  rcp<Epetra_Vector> vortWithGhost1_;
  rcp<Epetra_Vector> vortWithGhost2_;  
  // xgrid, ygrid
  rcp<Epetra_Vector> xx_;
  rcp<Epetra_Vector> yy_;
  // u,v fields 
  rcp<Epetra_Vector> UUNoGhost_;
  rcp<Epetra_Vector> VVNoGhost_;

  rcp<Epetra_Vector> streamNoGhost_;
  rcp<Epetra_Vector> streamWithGhost_;

  // poisson matrix since it remains the same  
  // what changes is only rhs
  rcp<Epetra_CrsMatrix> lapMat_;
};

  

}//end namespace apps

