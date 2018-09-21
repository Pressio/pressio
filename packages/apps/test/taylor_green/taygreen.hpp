
#ifndef TAYGREEN_HPP_
#define TAYGREEN_HPP_

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"
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
#include <map>


/* 2D TAYLOR GREEN problem, domain (0,1)^2
   mpi ranks are enumerated like this 

   * * *  
   0 1 2

   * * * 
   3 4 5 

   * * *
   6 7 8 

   indices: 
   i: is for rows, increases from top to bottom 
   j: is for cols, increases from left to right

   so rank 0 has indices (i,j) = (0,0)
   so rank 1 has indices (i,j) = (0,1)
   so rank 4 has indices (i,j) = (1,1)
 */


class TaylorGreen2d{
  
private:
  using nativeVec = Epetra_Vector;
  template<typename T>
  using rcp = std::shared_ptr<T>;

  using veci = std::vector<int>;
  using vecd = std::vector<double>;
  
public:
  using scalar_type = double;
  using state_type = Epetra_Vector;
  using space_residual_type = Epetra_Vector;

public:  
  TaylorGreen2d(int Ncell, Epetra_MpiComm * comm)
    : N_(Ncell), comm_(comm){}
  ~TaylorGreen2d() = default; 

  void setup(){
    storeRankInfo();
    storeRankNeighbors();
    printRankInfo();
    gridInfo();
    createMapNoGhost();
    createMapWithGhost();
    createImporters();
    createMesh();
    createFields();
    
    //create laplace operator
    lapMat_ = std::make_shared<Epetra_CrsMatrix>(Copy,*noGhMap_,5);
    fillLaplace();
    
    sten_.resize(4);
  };

public:
  void exchangeStreamFunctionGhosts(){
    psiWiGh_->Import(*psiNoGh_, *uniqueToGhostImporter_, Insert);
  }//end

  void exchangeVorticityGhosts(Epetra_Vector & omeUnique){
    omeWiGh_->Import(omeUnique, *uniqueToGhostImporter_, Insert);
  }//end

  int ijTok(int  i, int  j, int  N){
    return i*N + j;
  }
  

  void residual(Epetra_Vector & RHS,
		Epetra_Vector & ome,
		double t){
    
    solvePoisson( t, ome );
    exchangeStreamFunctionGhosts();
    computeVelocity();
    // // double errU = computeError(*UU_, t, "uu");
    // // double errV = computeError(*VV_, t, "vv");      
    exchangeVorticityGhosts(ome);

    double u,v,convX,convY,diffX,diffY;
    for (int i=0; i<lN_; i++){
     for (int j=0; j<lN_; j++){       
       int kk = ijTok(i, j, lN_);
       u = (*UU_)[kk];
       v = (*VV_)[kk];
       // compute stencil for vorticity
       computeStencil(*omeWiGh_, i, j);

       convX = u*( sten_[1] - sten_[0] )/(2.0*dx_); 
       convY = v*( sten_[2] - sten_[3] )/(2.0*dx_);
       diffX = ReInvDx_*( sten_[1] - 2*(*omeWiGh_)[kk] + sten_[0]);
       diffY = ReInvDx_*( sten_[2] - 2*(*omeWiGh_)[kk] + sten_[3]);
       RHS[kk] = -convX-convY + diffX + diffY;
     }
    }
    
  }//end
  

  void printField(std::string fileName,
		  const Epetra_Vector & f)
  {
    std::ofstream file;
    file.open( fileName );
    for(int i=0; i < lN_; i++){
      for(int j=0; j < lN_; j++){
       int kk = ijTok(i, j, lN_);
	file << std::fixed << std::setprecision(10) << f[kk] << " ";
      }
      file << std::endl;
    }
    file.close();    
  }
  //------------------------------------------------
  

  void run(double dt, int Nsteps = 1)
  {
    double t = 0.0;

    // init vorticity
    for (int i=0; i<lNsq_; i++)
      (*omeNoGh_)[i] = vort( (*xx_)[i], (*yy_)[i], t );

    // time loop
    for (int iStep=1; iStep<=Nsteps; iStep++)
    {
      *yn_ = *omeNoGh_;      
      //RKs1
      residual( *omeNoGhK1_, *yn_, t);
      for (int i=0; i<lNsq_; i++)
      	(*omeNoGh_)[i] = (*yn_)[i] + (2./3.) * dt * (*omeNoGhK1_)[i];
      t += (2./3.) * dt;

      //RKs2
      residual( *omeNoGhK2_, *omeNoGh_, t);
      for (int i=0; i<lNsq_; i++)
      	(*omeNoGh_)[i] = (*yn_)[i]
       + dt * ( 0.25*(*omeNoGhK1_)[i] + (3./4.)*(*omeNoGhK2_)[i] );
      t = static_cast<double>(iStep) * dt;
      
      // residual( *omeNoGhK1_, *omeNoGh_, t);
      // for (int i=0; i<lNsq_; i++)
      // 	(*omeNoGh_)[i] += dt * (*omeNoGhK1_)[i];
      // t = static_cast<double>(iStep) * dt;
      
      if (iStep % 10 == 0)
	  printField("C_" +
		     std::to_string(iStep) + ".txt", *omeNoGh_);
      if (myR_==0){
      	std::cout << " done with step: " << iStep
		  << " t = " << t << std::endl;
      }
      
      // //      double errPSI = computeError(*psiNoGh_, t, "sf");
      // double errVOR = computeError(*omeNoGh_, t, "vo");
      // if (myR_==0){
      // 	std::cout << " step: " << iStep
      // 		  << " errVO = " << std::setprecision(6)
      // 		  << std::sqrt(errVOR/(N_*N_))
      // 		  << std::endl;
      // }
    }

  }//end 

  

  void computeStencil(const Epetra_Vector & f,
		      int i, int j){
    // current ij -> local id
    int kk = ijTok(i, j, lN_);

    // left 
    sten_[0] = (j==0) ?
      f[lNsq_+i] //pick from ghost
      : f[kk-1]; // pick previous id

    //right 
    sten_[1] = (j==lN_-1) ?
      f[lNsq_+lN_+i]
      : f[kk+1];

    // top 
    sten_[2] = (i==0) ?
      f[lNsq_+2*lN_+j]
      : f[kk-lN_];

    // bottom
    sten_[3] = (i==lN_-1) ?
      f[lNsq_+3*lN_+j]
      : f[kk+lN_];
  }//end

  
  void computeVelocity(){
    for (int i=0; i<lN_; i++){
      for (int j=0; j<lN_; j++){
	int kk = ijTok(i, j, lN_);
	computeStencil(*psiWiGh_, i, j);
	(*UU_)[kk] = ( sten_[2] - sten_[3] )/(2.0*dx_);
	(*VV_)[kk] = (sten_[0] - sten_[1] )/(2.0*dx_);
     }
    }
  }//end
  
  
  void createImporters(){
    uniqueToGhostImporter_ =         /* targMap, sourceMap */
      std::make_shared<Epetra_Import>(*wiGhMap_, *noGhMap_);
    //uniqueToGhostImporter_->Print(std::cout);
  }//end 
  
  
private:

  // // TGV
  // double vort(double x, double y, double t)const {
  //   return -2.0*std::cos(x)*std::cos(y)*std::exp(-2.*t/Re_);}
  // double stream(double x, double y, double t)const {
  //   return -std::cos(x)*std::cos(y)*std::exp(-2.*t/Re_);}
  // double vel_U(double x, double y, double t){
  //   return std::sin(y)*std::cos(x)*std::exp(-2.*t/Re_);}
  // double vel_V(double x, double y, double t){    
  //   return -std::cos(y)*std::sin(x)*std::exp(-2.*t/Re_);}
  

  double vort(double x, double y, double t)const {
    double x2 = 2.*PI_*x ;
    double vortX = 0.05*2.*PI_*std::cos(x2);
    
    if (y <= 0.5){
      double y2 = 30.*(y-0.25);
      return vortX - 30./std::cosh(y2);
    }
    else {
      double y2 = 30.*(0.75-y);
      return vortX + 30./std::cosh(y2);
    }
  }
  
  double stream(double x, double y, double t)const {
    double x2 = 2.*PI_*x;
    double psix = (0.025/PI_) * std::cos(x2);
    
    double psiy = 0.0;
    if (y <= 0.5){
      double y2 = 30.*(y-0.25);
      psiy = (1./30.) * std::log( std::cosh(y2) );
    }
    else {
      double y2 = 30.*(0.75-y);
      psiy = -(1./30.) * std::log( std::cosh(y2) );
    }
    return psix+psiy;
  }

  double vel_U(double x, double y, double t){
    if (y <= 0.5){
      double y2 = 30.*(y-0.25);
      return std::tanh(y2);}
    else {
      double y2 = 30.*(0.75-y);
      return std::tanh(y2);}
  }

  double vel_V(double x, double y, double t){
    double x2 = 2.*PI_*x;
    return 0.05 * std::sin(x2);
  }


  double computeError(const Epetra_Vector & myF,
		       double time, std::string f){
    double err = 0.0, tmp = 0.0;
    for (int i=0; i<lNsq_; i++){
      if (f == "uu"){
	tmp = std::abs(myF[i]-vel_U((*xx_)[i], (*yy_)[i], time));
      }
      else if (f == "vv"){
	tmp = std::abs(myF[i]-vel_V((*xx_)[i], (*yy_)[i], time));
      }
      else if (f == "sf"){
	tmp = std::abs(myF[i]-stream((*xx_)[i], (*yy_)[i], time));
      }
      else if (f == "vo"){
	tmp = std::abs(myF[i]-vort((*xx_)[i], (*yy_)[i], time));
      }
      err += tmp*tmp;
      // if (myR_==0){
      // 	std::cout << std::setprecision(14)
      // 		  <<  (*xx_)[i] << " "
      // 		  <<  (*yy_)[i] << " "
      // 		  << (*UU_)[i] << " "
      // 		  << vel_U((*xx_)[i],(*yy_)[i],time)
      // 		  << std::endl;
      // }
    }

    double gErr = 0.0;
    comm_->SumAll(&err, &gErr, nProc_);
    return gErr;

  }//end 

  
  
  void solvePoisson(double time, Epetra_Vector omega)
  {
    // initialize stramfunction
    for (int i=0; i<lNsq_; i++){
      (*psiNoGh_)[i] = 0.; //stream((*xx_)[i],(*yy_)[i],0.0);
    }      

    // fix the value of psi at given point 
    // because equation is defined up to a constant
    if (myR_==0) {
      omega[0] = stream((*xx_)[0], (*yy_)[0], time);
      //(*omeNoGh_)[0] = stream((*xx_)[0], (*yy_)[0], time);
    }
    
    Epetra_LinearProblem Problem( lapMat_.get(),
				  psiNoGh_.get(),
				  &omega );
    AztecOO Solver(Problem);
    Solver.SetAztecOption( AZ_solver, AZ_GMRESR);
    Solver.SetAztecOption( AZ_conv, AZ_r0);
    //Solver.SetAztecOption( AZ_conv, AZ_noscaled);
    Solver.SetAztecOption( AZ_output, AZ_none );
    // Solver.SetAztecOption( AZ_precond, AZ_Jacobi );
    // Solver.SetAztecOption( AZ_poly_ord, 10 );

    // ================= MultiLevelOperator SECTION ========================
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
    ML_Epetra::MultiLevelOperator  MLPrec(ml_handle, *comm_, *noGhMap_, *noGhMap_);
    // ========== End of MultiLevelOperator SECTION ========================

    // //set this operator as preconditioner for AztecOO
    Solver.SetPrecOperator(&MLPrec);
    // solve
    Solver.Iterate(5000, 1e-13);
    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_handle);  
    //Solver.Iterate(5000, 1e-6);
    vecd status( AZ_STATUS_SIZE );
    Solver.GetAllAztecStatus( status.data() );
    double AztecOOExitStatus = status.operator[](1);
    assert( AztecOOExitStatus== AZ_normal );
    // std::cout << "AztecOOExitStatus = " << AztecOOExitStatus
    // 	      << " " << AZ_normal << std::endl;
  }//end
  


  void fillLaplace()
  {  
    veci Indices(4);
    vecd ValueDiag(1);
    veci IndicesDiag(1);
    double diagVal = 4.0/(dx_*dx_);
    double offdv = -1.0/(dx_*dx_);
    vecd Values{offdv, offdv, offdv, offdv};

    int westID, eastID, northID, southID, diagID;
    int localCount = 0;
    for (int li=0; li<lN_; li++){
      for (int lj=0; lj<lN_; lj++)
	{
	  // global ID of current point
	  int thisGID = myGel_[localCount]; 
	  int gi = grid_glob_id_to_globi(lN_, myR_,
					 pi_, pj_, thisGID);
	  int gj = grid_glob_id_to_globj(lN_, myR_,
					 pi_, pj_, thisGID);

	  if (thisGID != 0 ){
	  
	    westID = (gj==0) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::left),
				       li, lN_-1, lN_)
	      : (lj==0) ? thisGID-lNsq_+lN_-1 : thisGID-1;
	    eastID = (gj==N_-1) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::right),
				       li, 0, lN_)
	      : (lj==lN_-1) ? thisGID+lNsq_-lN_+1 : thisGID+1;

	    northID = (gi==0) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::top),
				       lN_-1, lj, lN_)
	      : (li==0) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::top),
				       lN_-1, lj, lN_)
	      : thisGID-lN_;

	    southID = (gi==N_-1) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::bottom),
				       0, lj, lN_)
	      : (li==lN_-1) ?
	    grid_loc_ij_to_glob_id(neighbors_.at(neighTag::bottom),
				       0, lj, lN_)
	      : thisGID+lN_;
	  
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
  

  
  void createFields(){
    
    omeNoGh_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    omeWiGh_  = std::make_shared<Epetra_Vector>(*wiGhMap_);
    omeNoGhK1_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    omeNoGhK2_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    yn_ = std::make_shared<Epetra_Vector>(*noGhMap_);

    psiNoGh_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    psiWiGh_  = std::make_shared<Epetra_Vector>(*wiGhMap_);
    UU_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    VV_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    
  }//end
  
  
  int grid_glob_id_to_globi(int lN, int procRank,
			    int procRanki,
			    int procRankj, int k){
    int localLinIndex = k-procRank*lN*lN;
    return localLinIndex/lN + procRanki*lN;
  }
  int grid_glob_id_to_globj(int lN, int procRank,
			    int procRanki,
			    int procRankj, int k){
    int localLinIndex = k-procRank*lN*lN;
    return localLinIndex % lN + procRankj*lN;
  }

  void createMesh(){
    xx_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    yy_ = std::make_shared<Epetra_Vector>(*noGhMap_);
    int k = 0;
    for (auto const & it : myGel_){
      double x = dx_*0.5 +
	dx_ * grid_glob_id_to_globj(lN_, myR_,
				    pi_, pj_, it);
      double y = L_ - dx_*0.5 -
	dx_ * grid_glob_id_to_globi(lN_, myR_,
				    pi_, pj_, it);
      (*xx_)[k] = x;
      (*yy_)[k] = y;
      k++;
    }
    //xx_->Print(std::cout);
    
  }//end

  
  void gridInfo(){
    // grid info
    dx_ = L_/static_cast<int>(N_);
    lN_ = static_cast<int>(N_)/static_cast<int>(Npi_);
    lNsq_ = lN_*lN_;
    totN_ = N_*N_;
    ReInvDx_ = 1.0/(Re_ * dx_ * dx_);
    
  }//end

  
  int grid_loc_ij_to_glob_id(int procRank, int li,
			       int lj, int lN){
    return procRank*lN*lN + li*lN + lj;
  }
  

  void createMapNoGhost(){ 
    NumMyElem_ = 0;
    for(int i1=0; i1<lN_; ++i1){
      for(int j1=0; j1<lN_; ++j1){
  	int gID = grid_loc_ij_to_glob_id(myR_, i1, j1, lN_);
  	myGel_.emplace_back(gID);
  	NumMyElem_++;
      }
    }
    noGhMap_ = std::make_shared<Epetra_Map>(-1, (int) NumMyElem_,
  					  myGel_.data(),
  					  0, *comm_ ) ;
    //noGhMap_->Print(std::cout);
    
  }//end

  void createMapWithGhost(){ 
    NumMyElemwGh_ = 0;
    for(int i1=0; i1<lN_; ++i1){
      for(int j1=0; j1<lN_; ++j1){
  	int gID = grid_loc_ij_to_glob_id(myR_, i1, j1, lN_);
  	myGelwGh_.emplace_back(gID);
  	NumMyElemwGh_++;
      }
    }

    //-------------------
    //store left ghosts
    //-------------------
    {
      // find which rank my neighbor has
      auto lr = neighbors_[left];
      // the glob id of the top-right corner of my left neighbor
      auto jTRL = grid_loc_ij_to_glob_id(lr, 0, lN_-1, lN_);
      for (int k=0; k<lN_; k++){
	myGelwGh_.emplace_back(jTRL+k*lN_);
	NumMyElemwGh_++;
    }}

    //-------------------
    //store right ghosts
    //-------------------
    {
      // find which rank my neighbor has
      auto lr = neighbors_[right];
      // the glob id of the top-left corner of my right neighbor
      auto jTRL = grid_loc_ij_to_glob_id(lr, 0, 0, lN_);
      for (int k=0; k<lN_; k++){
	myGelwGh_.emplace_back(jTRL+k*lN_);
	NumMyElemwGh_++;
    }}
    
    //-------------------
    //store top ghosts
    //-------------------
    {
      // find which rank my neighbor has
      auto lr = neighbors_[top];
      // the glob id of the bottom-left corner of my top neighbor
      auto jTRL = grid_loc_ij_to_glob_id(lr, lN_-1, 0, lN_);
      for (int k=0; k<lN_; k++){
	myGelwGh_.emplace_back(jTRL+k);
	NumMyElemwGh_++;
    }}

    //-------------------
    //store bottom ghosts
    //-------------------
    {
      // find which rank my neighbor has
      auto lr = neighbors_[bottom];
      // the glob id of the top-left corner of my bottom neighbor
      auto jTRL = grid_loc_ij_to_glob_id(lr, 0, 0, lN_);
      for (int k=0; k<lN_; k++){
	myGelwGh_.emplace_back(jTRL+k);
	NumMyElemwGh_++;
    }}
    

    wiGhMap_ = std::make_shared<Epetra_Map>(-1,
					    (int) NumMyElemwGh_,
					    myGelwGh_.data(),
					    0, *comm_ ) ;
    //wiGhMap_->Print(std::cout);
    
  }//end


  
  // N = tot number , k = enumerates proc, 0,1,2,3,...
  inline int rank_to_pi(int N, int k){ return k/N; }
  inline int rank_to_pj(int N, int k){ return k % N; }
  void rank_to_pij(int & i, int & j, int N, int k){
    i = rank_to_pi(N,k);  j = rank_to_pj(N,k);
  }//end 
 
  int pij_to_rank(int iin, int jin, int N){
    return iin*N + jin;
  }//end
  
  void printRankInfo(){
    auto leftRank = neighbors_[neighTag::left];
    auto rightRank = neighbors_[neighTag::right];
    auto topRank = neighbors_[neighTag::top];
    auto bottomRank = neighbors_[neighTag::bottom];
    
    std::cout << " r = " << myR_
	      << " nPi_ = " << Npi_ 
	      << " nPj_ = " << Npj_ 
	      << " myIDi_= " << pi_
	      << " myIDj_= " << pj_
	      << " l,r,t,b = "
	      << leftRank << " " << rightRank << " "
	      << topRank << " " << bottomRank
	      // << " lN_ = " << lN_
	      // << " dx_ = " << dx_
	      << std::endl;
  }//end 

  
  void storeRankInfo(){
    myR_ = comm_->MyPID();
    nProc_ = comm_->NumProc();
    Npi_ = std::sqrt<int>(nProc_);
    Npj_ = Npi_;
    rank_to_pij(pi_, pj_, Npj_, myR_);

  }//end 

  
  void storeRankNeighbors(){
    int leftRank = pj_==0 ?
      pij_to_rank(pi_, Npj_-1, Npj_) : myR_-1;
    int rightRank = pj_==Npj_-1 ?
      pij_to_rank(pi_, 0, Npj_) : myR_+1;
    int topRank = pi_==0 ?
      pij_to_rank(Npi_-1, pj_, Npj_) : myR_-Npj_;
    int bottomRank = pi_==Npi_-1 ?
      pij_to_rank(0, pj_, Npj_) : myR_+Npj_;
    
    neighbors_.insert( {neighTag::left,  leftRank} );
    neighbors_.insert( {neighTag::right, rightRank} );
    neighbors_.insert( {neighTag::top,   topRank} );
    neighbors_.insert( {neighTag::bottom,bottomRank} );
    
  }//end 

  
private:
  const double PI_ = 3.14159265358979323846;
  
  const scalar_type xL_ = 0.0; //left side of domain 
const scalar_type xR_ = 1.;//2.*PI_; // right side of domain
  scalar_type L_ = xR_-xL_;
  int N_; // # of cells along each axis
  scalar_type dx_; // cell size
  int lN_, lNsq_, totN_;
  // scalar_type dxInv_; // inv of cell size
  // const int nonZrPerRow_ = 2;

  Epetra_MpiComm * comm_;
  int myR_;
  enum neighTag {left,right,top,bottom};
  std::map<neighTag,int> neighbors_;
  
  // pi_, pj_ : identify a rank in the 2d rank arrangement
  // 	   0      pj_     Npj_
  // 	   -------------------
  //  0   |       |         |
  //      |   0   |    1    |
  // 	  |       |         |
  // pi_  -------------------
  //      |       |         |
  // 	  |   2   |    3    |
  // Npi_ |       |         |
  // 	   -------------------
  int pi_, pj_;
  int Npi_, Npj_; // tot proc along each axis
  int nProc_;
  
  rcp<Epetra_Map> noGhMap_;
  int NumMyElem_;
  std::vector<int> myGel_;
  
  rcp<Epetra_Map> wiGhMap_;
  int NumMyElemwGh_;
  std::vector<int> myGelwGh_;

  // mesh
  rcp<Epetra_Vector> xx_;
  rcp<Epetra_Vector> yy_;

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
  
  rcp<Epetra_Import> uniqueToGhostImporter_;
  //  rcp<Epetra_Export> ghostToUniqueExporter_;
  rcp<Epetra_CrsMatrix> lapMat_;

  vecd sten_;  
  const double Re_ = 10000.;
  const double ReInv_ = 1.0/Re_;
  double ReInvDx_;
  
};

#endif 
