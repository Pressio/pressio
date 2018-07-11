
#ifndef APPS_TEST_BURGERS1D_EIGEN_OBSERVER_HPP_
#define APPS_TEST_BURGERS1D_EIGEN_OBSERVER_HPP_

#include "apps_burgers1d_eigen.hpp"

struct snapshot_collector
{
private:
  using state_t = apps::burgers1dEigen::state_type;
  using matrix = std::vector<state_t>;  

  matrix snapshots_;
  size_t count_;

public:
  void operator()(size_t step,
		  double t,
		  const state_t & x){
    if (step % 50 ==0 ){
      snapshots_.emplace_back(x);
      count_++;
    }
  }

  size_t getCount(){
    return count_;
  };

  void printLastStoreState(){
    auto const & finalSol = snapshots_.back();
    std::cout << "Last stored state " << std::endl;
    for (int i=0; i<finalSol.size(); ++i)
      std::cout << finalSol(i)  << " ";
    std::cout << std::endl;
    std::cout << std::endl;
  } 

  /* void printToFile(){ */
  /*   std::ofstream file; */
  /*   file.open( "out.txt" ); */
  /*   for (size_t step=0; step<count_; ++step){ */
  /*     auto const & sol = snapshots_[step]; */
  /*     for (int i=0; i< sol.size(); ++i) */
  /* 	file << std::fixed << std::setprecision(10) << sol(i) << " ";  */
  /*     file << std::endl; */
  /*   } */
  /*   file.close();     */
  /* }; */

};//end class

#endif
