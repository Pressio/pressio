
#ifndef ROM_TEST_BURGERS1D_EIGEN_GALPRO_OBSERVER_HPP_
#define ROM_TEST_BURGERS1D_EIGEN_GALPRO_OBSERVER_HPP_


struct snaps_collector
{
  using native_state_t = apps::Burgers1dEpetra::state_type;
  using ode_state_t = core::Vector<native_state_t>;
  using snap_t = core::MultiVector<Epetra_MultiVector>;
  size_t count_;
  snap_t shots_;
  double dt_;
  
public:
  snaps_collector(int nimages, const Epetra_Map & rowmap, double dt)
    : count_(0), shots_(rowmap, nimages), dt_(dt){}
  ~snaps_collector() = default;
  
public:
  void operator()(size_t step,
		  double t,
		  ode_state_t y)
  {
    if ( std::abs(t - 1.0) < 1e-8 ||
	 std::abs(t - 3.0) < 1e-8 ||
	 std::abs(t - 5.0) < 1e-8 ||
	 std::abs(t - 8.0) < 1e-8 ||
	 std::abs(t - 10.0) < 1e-8 ||
	 std::abs(t - 15.0) < 1e-8 ||
	 std::abs(t - 20.0) < 1e-8 ||
	 std::abs(t - 25.0) < 1e-8 ||
	 std::abs(t - 30.0) < 1e-8 ||
	 std::abs(t - 32.0) < 1e-8){
      *(*shots_.data())(count_) = *y.data();
      count_++;
    }
  }

  void printAll(){
    shots_.data()->Print(std::cout << std::setprecision(12));
    // for(int i=0; i<snapshots_.rows(); i++){
    //   for (int j=0; j<count_; ++j)
    // 	std::cout << std::setprecision(14) << snapshots_(i,j) << " ";
    //   std::cout << std::endl;
    // }
  }    

  snap_t getSnaps(){
    return shots_;
  };
  
  // size_t getCount(){
  //   return count_;
  // };

  // void printLastStoreState(){
  //   // auto const & finalSol = snapshots_.back();
  //   // std::cout << "Last stored state " << std::endl;
  //   // for (decltype(finalSol.size()) i=0; i<finalSol.size(); ++i)
  //   //   std::cout << finalSol[i]  << " ";
  //   // std::cout << std::endl;
  //   // std::cout << std::endl;
  // } 

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
