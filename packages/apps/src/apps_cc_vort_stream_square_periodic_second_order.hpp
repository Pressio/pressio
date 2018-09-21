
#ifndef APP_CC_VORT_STREAM_SQUARE_PERIODIC_SECOND_ORDER_HPP_
#define APP_CC_VORT_STREAM_SQUARE_PERIODIC_SECOND_ORDER_HPP_

#include <map>
#include <memory>
#include <vector>
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"
#include "apps_square_periodic_cc_mesh.hpp"

namespace apps{ 
namespace impl{ 


class VortStCCSqPeriodicSecOrd{
private:
  template<typename T>
  using rcp = std::shared_ptr<T>;
  
 public:
  VortStCCSqPeriodicSecOrd(Epetra_MpiComm * comm, double Re)
    : mpiInfo_(comm), meshInfo_(&mpiInfo_), Re_(Re){}
  
  ~VortStCCSqPeriodicSecOrd() = default;
  //-----------------------------------

  void setup(){
    omeWiGh_  = std::make_shared<Epetra_Vector>(*meshInfo_.wiGhMap_);
    omeNoGh_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    omeNoGhK1_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    omeNoGhK2_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    yn_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);

    psiNoGh_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    psiWiGh_  = std::make_shared<Epetra_Vector>(*meshInfo_.wiGhMap_);
    UU_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    VV_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);

    b_ = std::make_shared<Epetra_Vector>(*meshInfo_.noGhMap_);
    
    //create laplace operator
    lapMat_ = std::make_shared<Epetra_CrsMatrix>(Copy,*meshInfo_.noGhMap_,5);
    fillLaplace();
    
    ReInv_ = 1.0/Re_;
    ReInvDx_ = 1.0/(Re_ * meshInfo_.dx_ * meshInfo_.dx_);
  }//
  //-----------------------------------


  virtual void residual(Epetra_Vector & RHS, Epetra_Vector & ome, double t){

    // rhs is the current vorticity
    *b_ = ome;
    // initialize stramfunction
    for (int i=0; i<meshInfo_.lNsq_; i++){
      (*psiNoGh_)[i] = 0.0; //stream((*xx_)[i],(*yy_)[i],0.0);
    }      
    // fix the value of RHS at given point 
    // because equation is defined up to a constant
    fixRHSForPoissonSolve(t);
    
    solvePoisson();
    exchangeStreamFunctionGhosts();
    computeVelocity(t);
    exchangeVorticityGhosts(ome);

    double u,v,convX,convY,diffX,diffY;
    for (int i=0; i<meshInfo_.lN_; i++){
     for (int j=0; j<meshInfo_.lN_; j++){
       int kk = meshInfo_.local_ijTok(i, j);
       u = (*UU_)[kk];
       v = (*VV_)[kk];
       // compute stencil for vorticity
       computeStencil(*omeWiGh_, i, j);

       convX = u*( sten_[1] - sten_[0] )/(2.0*meshInfo_.dx_); 
       convY = v*( sten_[2] - sten_[3] )/(2.0*meshInfo_.dx_);
       diffX = ReInvDx_*( sten_[1] - 2*(*omeWiGh_)[kk] + sten_[0]);
       diffY = ReInvDx_*( sten_[2] - 2*(*omeWiGh_)[kk] + sten_[3]);
       RHS[kk] = -convX-convY + diffX + diffY;
     }
    }
  }//end
  //-----------------------------------

  virtual void fixRHSForPoissonSolve(double time) = 0;
  //-----------------------------------

  virtual void initializeVorticity() = 0;
  //-----------------------------------
  
  void computeStencil(const Epetra_Vector & F, int i, int j){
    meshInfo_.stencilIDs(i,j, lid, rid, tid, bid);
    // left 
    sten_[0] = F[lid];
    //right 
    sten_[1] = F[rid];
    // top 
    sten_[2] = F[tid];
    // bottom
    sten_[3] = F[bid];
  }//end
  //-----------------------------------


  void printField(std::string fileName,
		  const Epetra_Vector & f){
    std::ofstream file;
    file.open( fileName );
    for(int i=0; i < meshInfo_.lN_; i++){
      for(int j=0; j < meshInfo_.lN_; j++){
       int kk = meshInfo_.local_ijTok(i, j);
	file << std::fixed
	     << std::setprecision(10) << f[kk] << " ";
      }
      file << std::endl;
    }
    file.close();    
  }
  //------------------------------------------------
  
  
  virtual void computeVelocity(double time){
    for (int i=0; i<meshInfo_.lN_; i++){
      for (int j=0; j<meshInfo_.lN_; j++){
  	int kk = meshInfo_.local_ijTok(i, j);
  	computeStencil(*psiWiGh_, i, j);
  	(*UU_)[kk] = ( sten_[2] - sten_[3] )/(2.0*meshInfo_.dx_);
  	(*VV_)[kk] = (sten_[0] - sten_[1] )/(2.0*meshInfo_.dx_);
     }
    }
  }//end
  //-----------------------------------

  
  void exchangeStreamFunctionGhosts(){
    psiWiGh_->Import(*psiNoGh_,*meshInfo_.uniqueToGhostImporter_,Insert);
  }//end
  //-----------------------------------

  void exchangeVorticityGhosts(Epetra_Vector & omeUnique){
    omeWiGh_->Import(omeUnique,*meshInfo_.uniqueToGhostImporter_,Insert);
  }//end
  //-----------------------------------

  
  void solvePoisson()
  {
    Epetra_LinearProblem Problem(lapMat_.get(), psiNoGh_.get(), b_.get() );
    AztecOO Solver(Problem);
    Solver.SetAztecOption( AZ_solver, AZ_GMRESR);
    Solver.SetAztecOption( AZ_conv, AZ_r0);
    //Solver.SetAztecOption( AZ_conv, AZ_noscaled);
    Solver.SetAztecOption( AZ_output, AZ_none );
    // Solver.SetAztecOption( AZ_precond, AZ_Jacobi );
    // Solver.SetAztecOption( AZ_poly_ord, 10 );

    // ========== MultiLevelOperator SECTION ========================
    int nLevels = 10;            // maximum number of levels
    int maxMgLevels = 6;         //
    ML_Set_PrintLevel(0);       // print level (0 silent, 10 verbose)
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
    nLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle,
  						  maxMgLevels-1,
  						  ML_DECREASING,
  						  agg_object);
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
    ML_Epetra::MultiLevelOperator  MLPrec(ml_handle, *mpiInfo_.comm_,
  					  *meshInfo_.noGhMap_,
  					  *meshInfo_.noGhMap_);
    // ====== End of MultiLevelOperator SECTION ===============

    // //set this operator as preconditioner for AztecOO
    Solver.SetPrecOperator(&MLPrec);
    // solve
    Solver.Iterate(5000, 1e-13);
    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_handle);  
    //Solver.Iterate(5000, 1e-6);
    std::vector<double> status( AZ_STATUS_SIZE );
    Solver.GetAllAztecStatus( status.data() );
    double AztecOOExitStatus = status.operator[](1);
    assert( AztecOOExitStatus== AZ_normal );
    // std::cout << "AztecOOExitStatus = " << AztecOOExitStatus
    // 	      << " " << AZ_normal << std::endl;
  }//end
  //-----------------------------------


  
  void fillLaplace()
  {
    double dx = meshInfo_.dx_;
    int lN = meshInfo_.lN_;
    int lNsq = lN*lN;
    auto myGel = meshInfo_.myGel_;
    
    using veci = std::vector<int>;
    using vecd = std::vector<double>;
    
    veci Indices(4);
    vecd ValueDiag(1);
    veci IndicesDiag(1);
    double diagVal = 4.0/(dx*dx);
    double offdv = -1.0/(dx*dx);
    vecd Values{offdv, offdv, offdv, offdv};

    int westID, eastID, northID, southID;//, diagID;
    int localCount = 0;
    for (int li=0; li<lN; li++){
      for (int lj=0; lj<lN; lj++)
  	{
  	  // global ID of current point
  	  int thisGID = myGel[localCount]; 
  	  int gi = meshInfo_.glob_id_to_globi(thisGID);
  	  int gj = meshInfo_.glob_id_to_globj(thisGID);

  	  if (thisGID != 0 ){
	  
      westID = (gj==0) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.leftNeighbor(), li, lN-1)
  	: (lj==0) ? thisGID-lNsq+lN-1 : thisGID-1;
      eastID = (gj==meshInfo_.N_-1) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.rightNeighbor(),li, 0)
  	: (lj==lN-1) ? thisGID+lNsq-lN+1 : thisGID+1;

      northID = (gi==0) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.topNeighbor(),lN-1, lj)
  	: (li==0) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.topNeighbor(),lN-1, lj)
  	: thisGID-lN;

      southID = (gi==meshInfo_.N_-1) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.bottomNeighbor(),0, lj)
  	: (li==lN-1) ?
  	meshInfo_.loc_ij_to_glob_id( mpiInfo_.bottomNeighbor(),0, lj)
  	: thisGID+lN;

  	    Indices[0] = westID; 
  	    Indices[1] = eastID;
  	    Indices[2] = northID;
  	    Indices[3] = southID;
  	    //off diagonal
  	    lapMat_->InsertGlobalValues(thisGID, 4,
  				 Values.data(), Indices.data() );
  	    // add diagonal
  	    ValueDiag[0] = diagVal;
  	    IndicesDiag[0] = thisGID;
  	    lapMat_->InsertGlobalValues( thisGID, 1,
  				  ValueDiag.data(),
  				  IndicesDiag.data() );
	    
  	  }
  	  else
  	    {
  	      ValueDiag[0] = 1.0;
  	      IndicesDiag[0] = thisGID;
  	      lapMat_->InsertGlobalValues( thisGID, 1,
  				    ValueDiag.data(),
  				    IndicesDiag.data() );
  	    }
	
  	  localCount++;
  	}//end j
    }//end i

    lapMat_->FillComplete();
    //lapMat_->Print(std::cout);
  }//edn method 
  //-----------------------------------

  
protected:
  apps::impl::MpiSquareGrid mpiInfo_;
  apps::impl::SquarePeriodicCCMesh meshInfo_;
  
  double Re_;
  double ReInv_;
  double ReInvDx_;

  // vorticity
  rcp<Epetra_Vector> yn_;
  rcp<Epetra_Vector> omeNoGh_;
  rcp<Epetra_Vector> omeWiGh_;
  rcp<Epetra_Vector> omeNoGhK1_;
  rcp<Epetra_Vector> omeNoGhK2_;
  // psi
  rcp<Epetra_Vector> psiNoGh_;
  rcp<Epetra_Vector> psiWiGh_;
  // u,v 
  rcp<Epetra_Vector> UU_;
  rcp<Epetra_Vector> VV_;

  rcp<Epetra_Vector> b_;//rhs for poisson
  
  rcp<Epetra_CrsMatrix> lapMat_;

  int lid, rid, tid, bid;
  std::array<double,4> sten_;
  
};// end 

} //end namespace impl
} //end namespace apps
#endif



  // void residual(Epetra_Vector & RHS,
  // 		Epetra_Vector & ome,
  // 		double t){

  //   // fill psi and omega
  //   PoissonSolver::solve( *lapMat_ );
  //   exchangeStreamFunctionGhosts();
  //   computeVelocity();
  //   // // double errU = computeError(*UU_, t, "uu");
  //   // // double errV = computeError(*VV_, t, "vv");      
  //   exchangeVorticityGhosts(ome);

  //   // double u,v,convX,convY,diffX,diffY;
  //   // for (int i=0; i<lN; i++){
  //   //  for (int j=0; j<lN; j++){       
  //   //    int kk = ijTok(i, j, lN);
  //   //    u = (*UU_)[kk];
  //   //    v = (*VV_)[kk];
  //   //    // compute stencil for vorticity
  //   //    computeStencil(*omeWiGh_, i, j);
  //   //    convX = u*( sten_[1] - sten_[0] )/(2.0*dx_); 
  //   //    convY = v*( sten_[2] - sten_[3] )/(2.0*dx_);
  //   //    diffX = ReInvDx_*( sten_[1] - 2*(*omeWiGh_)[kk] + sten_[0]);
  //   //    diffY = ReInvDx_*( sten_[2] - 2*(*omeWiGh_)[kk] + sten_[3]);
  //   //    RHS[kk] = -convX-convY + diffX + diffY;
  //   //  }
  //   // } 
  // }//end
  // //-----------------------------------
