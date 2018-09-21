
#ifndef APP_VORT_STREAM_AUX_HPP_
#define APP_VORT_STREAM_AUX_HPP_

#include "Epetra_MpiComm.h"
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"
#include <map>
#include <memory>
#include <vector>
#include "Epetra_Vector.h"


namespace apps{ 
namespace impl{ 


class Mesh
{
private:
  template<typename T>
  using rcp = std::shared_ptr<T>;

public:
  Mesh(apps::impl::MpiSquareGrid * mpiInfo)
    : mpiInfo_(mpiInfo){}
  ~Mesh() = default;

  void setup(double xL, double xR, int N){
    myR_ = mpiInfo_->myR_;
    N_ = N;
    xL_ = xL;
    xR_ = xR;    
    L_ = xR_ - xL_;
    dx_ = L_/static_cast<int>(N_);
    dxInv_ = 1.0/dx_;
    lN_ = static_cast<int>(N_)/static_cast<int>(mpiInfo_->Npi_);
    lNsq_ = lN_*lN_;
    totN_ = N_*N_;

    createMapNoGhost();
    createMapWithGhost();
    createMesh();
  }


  void fillLaplace(Epetra_CrsMatrix & A, bool invertSign = false)
  {
    using veci = std::vector<int>;
    using vecd = std::vector<double>;
    
    veci Indices(4);
    vecd ValueDiag(1);
    veci IndicesDiag(1);
    double diagVal = -4.0/(dx_*dx_);
    double offdv = 1.0/(dx_*dx_);
    if (invertSign){
      diagVal *= -1.0;
      offdv *= -1.0;
    }
    vecd Values{offdv, offdv, offdv, offdv};

    int westID, eastID, northID, southID, diagID;
    int localCount = 0;
    for (int li=0; li<lN_; li++){
      for (int lj=0; lj<lN_; lj++)
	{
	  // global ID of current point
	  int thisGID = myGel_[localCount]; 
	  int gi = grid_glob_id_to_globi(lN_, myR_,
					 mpiInfo_->pi_,
					 mpiInfo_->pj_, thisGID);
	  int gj = grid_glob_id_to_globj(lN_, myR_,
					 mpiInfo_->pi_,
					 mpiInfo_->pj_, thisGID);

	  if (thisGID != 0 ){
	  
	    westID = (gj==0) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->leftNeighbor(),
				     li, lN_-1, lN_)
	      : (lj==0) ? thisGID-lNsq_+lN_-1 : thisGID-1;
	    eastID = (gj==N_-1) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->rightNeighbor(),
				     li, 0, lN_)
	      : (lj==lN_-1) ? thisGID+lNsq_-lN_+1 : thisGID+1;

	    northID = (gi==0) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->topNeighbor(),
				     lN_-1, lj, lN_)
	      : (li==0) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->topNeighbor(),
				     lN_-1, lj, lN_)
	      : thisGID-lN_;

	    southID = (gi==N_-1) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->bottomNeighbor(),
				     0, lj, lN_)
	      : (li==lN_-1) ?
	      grid_loc_ij_to_glob_id( mpiInfo_->bottomNeighbor(),
				     0, lj, lN_)
	      : thisGID+lN_;
	  
	    Indices[0] = westID; 
	    Indices[1] = eastID;
	    Indices[2] = northID;
	    Indices[3] = southID;
	    //off diagonal
	    A.InsertGlobalValues(thisGID, 4,
					Values.data(), Indices.data() );
	    // add diagonal
	    ValueDiag[0] = diagVal;
	    IndicesDiag[0] = thisGID;
	    A.InsertGlobalValues( thisGID, 1,
					 ValueDiag.data(),
					 IndicesDiag.data() );
	    
	  }
	  else
	    {
	      ValueDiag[0] = 1.0;
	      IndicesDiag[0] = thisGID;
	      A.InsertGlobalValues( thisGID, 1,
					   ValueDiag.data(),
					   IndicesDiag.data() );
	    }
	
	  localCount++;
	}//end j
    }//end i

    A.FillComplete();
    //lapMat_->Print(std::cout);
  }//edn method 


  void stencilIDs(int i, int j,
		  int & lid, int & rid, int & tid, int & bid)
  {
    // current ij -> local id
    int kk = local_ijTok(i, j, lN_);
    
    lid = (j==0) ? lNsq_+i //pick from ghost
      : kk-1; // pick previous id
    //right 
    rid = (j==lN_-1) ? lNsq_+lN_+i : kk+1;
    // top 
    tid = (i==0) ? lNsq_+2*lN_+j : kk-lN_;
    // bottom
    bid = (i==lN_-1) ? lNsq_+3*lN_+j : kk+lN_;    
  }
  
  
private:
  
  int local_ijTok(int  i, int  j, int  N){
    return i*N + j;
  }

  int grid_loc_ij_to_glob_id(int procRank, int li, int lj, int lN){
    return procRank*lN*lN + li*lN + lj;
  }

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
				    mpiInfo_->pi_,
				    mpiInfo_->pj_, it);
      double y = L_ - dx_*0.5 -
	dx_ * grid_glob_id_to_globi(lN_, myR_,
				    mpiInfo_->pi_,
				    mpiInfo_->pj_, it);
      (*xx_)[k] = x;
      (*yy_)[k] = y;
      k++;
    }
    //xx_->Print(std::cout);
  }//end

  
  void createImporters(){
    uniqueToGhostImporter_ =         /* targMap, sourceMap */
      std::make_shared<Epetra_Import>(*wiGhMap_, *noGhMap_);
    //uniqueToGhostImporter_->Print(std::cout);
  }//end 

  
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
					    0, *mpiInfo_->comm_ ) ;
    //noGhMap_->Print(std::cout);    
  }//end
  //--------------------------------------

  void createMapWithGhost()
  { 
    NumMyElemwGh_ = 0;
    for(int i1=0; i1<lN_; ++i1){
      for(int j1=0; j1<lN_; ++j1){
  	int gID = grid_loc_ij_to_glob_id(myR_, i1, j1, lN_);
  	myGelwGh_.emplace_back(gID);
  	NumMyElemwGh_++;
      }
    }
    storeLeftGhostInfo();
    storeRightGhostInfo();
    storeTopGhostInfo();
    storeBottomGhostInfo();

    wiGhMap_ = std::make_shared<Epetra_Map>(-1, (int) NumMyElemwGh_,
					    myGelwGh_.data(),
					    0, *mpiInfo_->comm_ ) ;
    //wiGhMap_->Print(std::cout); 
  }//end
  //---------------------------------------------

  void storeLeftGhostInfo(){
    //store left ghosts
    // find which rank my neighbor has
    auto lr = mpiInfo_->leftNeighbor();
    // the glob id of the top-right corner of my left neighbor
    auto jTRL = grid_loc_ij_to_glob_id(lr, 0, lN_-1, lN_);
    for (int k=0; k<lN_; k++){
      myGelwGh_.emplace_back(jTRL+k*lN_);
      NumMyElemwGh_++;
    }
  }//end

  void storeRightGhostInfo(){
    //store right ghosts
    // find which rank my neighbor has
    auto lr = mpiInfo_->rightNeighbor();
    // the glob id of the top-left corner of my right neighbor
    auto jTRL = grid_loc_ij_to_glob_id(lr, 0, 0, lN_);
    for (int k=0; k<lN_; k++){
      myGelwGh_.emplace_back(jTRL+k*lN_);
      NumMyElemwGh_++;
    }
  }//end

  void storeTopGhostInfo(){
    //store top ghosts
    // find which rank my neighbor has
    auto lr = mpiInfo_->topNeighbor();
    // the glob id of the bottom-left corner of my top neighbor
    auto jTRL = grid_loc_ij_to_glob_id(lr, lN_-1, 0, lN_);
    for (int k=0; k<lN_; k++){
      myGelwGh_.emplace_back(jTRL+k);
      NumMyElemwGh_++;
    }
  }//end

  void storeBottomGhostInfo(){
    //store bottom ghosts
    // find which rank my neighbor has
    auto lr = mpiInfo_->bottomNeighbor();
    // the glob id of the top-left corner of my bottom neighbor
    auto jTRL = grid_loc_ij_to_glob_id(lr, 0, 0, lN_);
    for (int k=0; k<lN_; k++){
      myGelwGh_.emplace_back(jTRL+k);
      NumMyElemwGh_++;
    }
  }//end
    


  
public:
  apps::impl::MpiSquareGrid * mpiInfo_;
  int myR_;

  double xL_;
  double xR_;
  double L_;
  
  int N_; // # of cells along each axis
  double dx_; // cell size
  int lN_, lNsq_, totN_;
  double dxInv_; // inv of cell size

  rcp<Epetra_Map> noGhMap_;
  int NumMyElem_;
  std::vector<int> myGel_;
  
  rcp<Epetra_Map> wiGhMap_;
  int NumMyElemwGh_;
  std::vector<int> myGelwGh_;

  rcp<Epetra_Import> uniqueToGhostImporter_;
  
  // mesh
  rcp<Epetra_Vector> xx_;
  rcp<Epetra_Vector> yy_;
  
};// end mesh
////////////////////////////////////





////////////////////////////////////
class Field{
private:
  template<typename T>
  using rcp = std::shared_ptr<T>;
  
 public:
  Field() = default;
  ~Field() = default;


  void storeGridInfo(apps::impl::Mesh & grid_){
    grid_ = &grid_;
  }
  
  void computeStencil(int i, int j, std::vector<double> & sten){
    grid_->stencilIDs(i,j, lid, rid, tid, bid);
    // left 
    sten_[0] = ff[lid];
    //right 
    sten_[1] = ff[rid];
    // top 
    sten_[2] = ff[tid];
    // bottom
    sten_[3] = ff[bid];
  }//end
  
  void import(Epetra_Vector & sourceF, Epetra_Import & importer){
    //psiWiGh_->Import(*psiNoGh_, *uniqueToGhostImporter_, Insert);
    ff_->Import(sourceF, importer, Insert);
  }//end

  
  void create(const Epetra_Map & map){
    ff_ = std::make_shared<Epetra_Vector>(map);
  }//end 
  
public:
  rcp<Epetra_Vector> ff_;
  apps::impl::Mesh * grid_;
  int lid, rid, tid, bid;
  
};// end fields
////////////////////////////////////
  




  
class PoissonSolver{
  

};




  // void printField(std::string fileName,
  // 		  const Epetra_Vector & f)
  // {
  //   std::ofstream file;
  //   file.open( fileName );
  //   for(int i=0; i < lN_; i++){
  //     for(int j=0; j < lN_; j++){
  //      int kk = ijTok(i, j, lN_);
  // 	file << std::fixed << std::setprecision(10) << f[kk] << " ";
  //     }
  //     file << std::endl;
  //   }
  //   file.close();    
  // }
  // //------------------------------------------------



  
  

} //end namespace impl
} //end namespace apps
#endif 
