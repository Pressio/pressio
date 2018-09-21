
#ifndef APPS_SQUARE_PERIODIC_CC_MESH_HPP__
#define APPS_SQUARE_PERIODIC_CC_MESH_HPP__

#include <map>
#include <memory>
#include <vector>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Vector.h>
#include "apps_mpi_square_grid.hpp"

namespace apps{ 
namespace impl{ 
  
class SquarePeriodicCCMesh{  
private:
  template<typename T>
  using rcp = std::shared_ptr<T>;

public:
  SquarePeriodicCCMesh(apps::impl::MpiSquareGrid * mpiInfo)
    : mpiInfo_(mpiInfo){}
  ~SquarePeriodicCCMesh() = default;

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
    createImporters();
    createMesh();
  }

  void stencilIDs(int i, int j,
		  int & lid, int & rid, int & tid, int & bid)
  {
    // current ij -> local id
    int kk = local_ijTok(i, j, lN_);

    // left
    lid = (j==0) ? lNsq_+i //pick from ghost
      : kk-1; // pick previous id
    //right 
    rid = (j==lN_-1) ? lNsq_+lN_+i : kk+1;
    // top 
    tid = (i==0) ? lNsq_+2*lN_+j : kk-lN_;
    // bottom
    bid = (i==lN_-1) ? lNsq_+3*lN_+j : kk+lN_;    
  }

  
  int glob_id_to_globi(int gID){
    int localLinIndex = gID-myR_*lN_*lN_;
    return localLinIndex/lN_ + mpiInfo_->pi_*lN_;
  }
  
  int glob_id_to_globj(int gID){
    int localLinIndex = gID-myR_*lN_*lN_;
    return localLinIndex % lN_ + mpiInfo_->pj_*lN_;
  }

  int loc_ij_to_glob_id(int procRank, int li, int lj){
    return procRank*lN_*lN_ + li*lN_ + lj;
  }

  int local_ijTok(int  i, int  j){
    return i*lN_ + j;
  }
  
private:
  
  int local_ijTok(int  i, int  j, int  N){
    return i*N + j;
  }

  int loc_ij_to_glob_id(int procRank, int li, int lj, int lN){
    return procRank*lN*lN + li*lN + lj;
  }

  int glob_id_to_globi(int lN, int procRank,
			    int procRanki,
			    int procRankj, int k){
    int localLinIndex = k-procRank*lN*lN;
    return localLinIndex/lN + procRanki*lN;
  }
  
  int glob_id_to_globj(int lN, int procRank,
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
	dx_ * glob_id_to_globj(lN_, myR_,
			       mpiInfo_->pi_,
			       mpiInfo_->pj_, it);
      double y = L_ - dx_*0.5 -
	dx_ * glob_id_to_globi(lN_, myR_,
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
  	int gID = loc_ij_to_glob_id(myR_, i1, j1, lN_);
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
  	int gID = loc_ij_to_glob_id(myR_, i1, j1, lN_);
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
    auto jTRL = loc_ij_to_glob_id(lr, 0, lN_-1, lN_);
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
    auto jTRL = loc_ij_to_glob_id(lr, 0, 0, lN_);
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
    auto jTRL = loc_ij_to_glob_id(lr, lN_-1, 0, lN_);
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
    auto jTRL = loc_ij_to_glob_id(lr, 0, 0, lN_);
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

} //end namespace impl
} //end namespace apps
#endif
