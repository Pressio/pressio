
#ifndef APP_MPI_SQUARE_GRID_HPP_
#define APP_MPI_SQUARE_GRID_HPP_

#include "Epetra_MpiComm.h"
#include <map>
#include <memory>
#include <vector>

namespace apps{ 
namespace impl{ 

  
class MpiSquareGrid{
 public:
  MpiSquareGrid(Epetra_MpiComm * comm) : comm_(comm){
    storeRankInfo();
    storeRankNeighbors();
    printRankInfo();    
  }
  ~MpiSquareGrid() = default;

  int leftNeighbor(){
    return neighbors_[left];}
  int rightNeighbor(){
    return neighbors_[right];}
  int topNeighbor(){
    return neighbors_[top];}
  int bottomNeighbor(){
    return neighbors_[bottom];}
  
private:
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

private:
  inline int rank_to_pi(int N, int k){ 
    return k/N; 
  }//end
  
  inline int rank_to_pj(int N, int k){ 
    return k % N; 
  }//end 

  inline void rank_to_pij(int & i, int & j, int N, int k){
    i = rank_to_pi(N,k);  
    j = rank_to_pj(N,k);
  }//end 
 
  inline int pij_to_rank(int iin, int jin, int N){
    return iin*N + jin;
  }//end

public:
  Epetra_MpiComm * comm_;
  int myR_; // my rank
  enum neighTag {left,right,top,bottom};
  std::map<neighTag,int> neighbors_;
  
  // pi_, pj_ : identify a rank in the 2d rank arrangement
  //     0      pj_     Npj_
  //     -------------------
  //  0   |       |         |
  //      |   0   |    1    |
  //      |       |         |
  // pi_  -------------------
  //      |       |         |
  //      |   2   |    3    |
  // Npi_ |       |         |
  //     -------------------
  int pi_, pj_;
  int Npi_, Npj_; // tot proc along each axis
  int nProc_;
  
};// end MpiSquareGrid

} //end namespace impl
} //end namespace apps
#endif
