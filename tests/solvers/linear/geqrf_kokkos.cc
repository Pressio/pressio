
#include "pressio_solvers.hpp"

int main(int argc, char *argv[])
{

  // device layout
  using d_layout = Kokkos::LayoutLeft;
  // // host layout
  // using h_layout = Kokkos::LayoutLeft;

  using exe_space = Kokkos::DefaultExecutionSpace;

  // typedef for vector
  using k1d_d = Kokkos::View<double*, d_layout, exe_space>;
  using k1d_h = typename k1d_d::host_mirror_type;

  using wv = ::pressio::containers::Vector<k1d_d>;
  static_assert( ::pressio::containers::details::traits<wv>::is_static == 0, "" );
  static_assert( ::pressio::containers::predicates::is_vector_wrapper_kokkos<wv>::value, "" );

  // typedef for matrix
  using k2d_d = Kokkos::View<double**, d_layout, exe_space>;
  using k2d_h = typename k2d_d::host_mirror_type;

  using wmat = ::pressio::containers::DenseMatrix<k2d_d>;
  static_assert( ::pressio::containers::details::traits<wmat>::is_static == 0, "" );
  static_assert( ::pressio::containers::predicates::is_dense_matrix_wrapper_kokkos<wmat>::value, "" );

  // size of problems, rows and cols
  constexpr int Nr = 4;
  constexpr int Nc = 4;

  Kokkos::initialize (argc, argv);
  {
  // ----------------------
  // create/fill the matrix
  // ----------------------
  // create host matrix
  k2d_h A_h("Ah", Nr, Nc);
  // fill host matrix
  A_h(0,0)=1.; A_h(0,1)=3.; A_h(0,2)=2.; A_h(0,3)=3.;
  A_h(1,0)=2.; A_h(1,1)=4.; A_h(1,2)=1.; A_h(1,3)=2.;
  A_h(2,0)=1.; A_h(2,1)=6.; A_h(2,2)=3.; A_h(2,3)=1.;
  A_h(3,0)=0.; A_h(3,1)=2.; A_h(3,2)=2.; A_h(3,3)=1.;

  // create device matrix
  k2d_d A_d("Ad", Nr, Nc);
  // deep copy host -> device
  Kokkos::deep_copy(A_d, A_h);

  // ------------------
  // Construct rhs
  // ------------------
  // create host
  k1d_h b_h("bh", Nr);
  // fill host matrix
  b_h(0) = 4.;
  b_h(1) = 1.;
  b_h(2) = 3.0;
  b_h(3) = 2.;

  // create device
  k1d_d b_d("bd", Nr);
  // copy host -> device
  Kokkos::deep_copy(b_d, b_h);

  //--------------------
  // solve system
  //--------------------
  // wrap device matrix and vectors
  wmat wA_d(A_d);
  wv wb_d(b_d);
  wv wx_d("x", Nc);

  // linear solver
  using solver_tag   = pressio::solvers::linear::direct::geqrf;
  using linear_solver_t = pressio::solvers::linear::Solver<solver_tag, wmat>;
  linear_solver_t lsObj;
  lsObj.solveAllowMatOverwrite(wA_d, wb_d, wx_d);

  // create host of solution
  k1d_h x_h("xh", Nc);
  Kokkos::deep_copy(x_h, *wx_d.data());

  std::cout << x_h(0) << " "
  	    << x_h(1) << " "
  	    << x_h(2) << " "
  	    << x_h(3) << std::endl;

  std::vector<double> sum(Nc, 0.0);
  for (int i=0; i<Nr; ++i)
    for (int j=0; j<Nc; ++j)
      sum[i] += A_h(i,j) * x_h(j);

  std::string strOutput  = "PASSED";
  for (int i=0; i<Nr; ++i){
    const auto err = sum[i] - b_h(i);
    if (err > 1e-12) strOutput = "FAILED";
    // EXPECT_NEAR( b_h(i), sum[i], 1e-12);
  }
  std::cout << strOutput << std::endl; 
  }

  Kokkos::finalize();
  return 0;
}
