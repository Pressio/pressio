
#ifndef ROM_TEST_BURGERS1D_EIGEN_GALPRO_OBSERVER_HPP_
#define ROM_TEST_BURGERS1D_EIGEN_GALPRO_OBSERVER_HPP_

struct snapshot_collector
{
  using native_state_t = apps::burgers1dEigen::state_type;
  using state_t = core::vector<native_state_t>;
  using matrix = Eigen::MatrixXd;
  matrix snapshots_;
  size_t count_;

public:
  snapshot_collector(int ncell, int nimages) : count_(0){
    snapshots_.resize(ncell, nimages);
  }
  
public:
  void operator()(size_t step,
		  double t,
		  state_t y){
    // //    if (step % 50 ==0 || step==0){
    //   for(int i=0; i<y.size(); i++){
    // 	snapshots_(i,count_) = y[i];
    //   }
    //   count_++;
    //   //}
  }

  void printAll()
  { 
    for(int i=0; i<snapshots_.rows(); i++){
      for (int j=0; j<count_; ++j)
	std::cout << std::setprecision(14) << snapshots_(i,j) << " ";
      std::cout << std::endl;
    }
  }    

  size_t getCount(){
    return count_;
  };

  void printLastStoreState(){
    // auto const & finalSol = snapshots_.back();
    // std::cout << "Last stored state " << std::endl;
    // for (decltype(finalSol.size()) i=0; i<finalSol.size(); ++i)
    //   std::cout << finalSol[i]  << " ";
    // std::cout << std::endl;
    // std::cout << std::endl;
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
