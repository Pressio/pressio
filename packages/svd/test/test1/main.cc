
#include "svd_solver_eigen.hpp"
#include "svd_solver_traits.hpp"

#include "vector/core_vector_epetra.hpp"
#include "vector/core_vector_serial_arbitrary.hpp"
#include "vector/core_vector_std.hpp"
#include "vector/core_vector_eigen.hpp"

#include "matrix/core_matrix_traits.hpp"
#include "matrix/core_matrix_eigen.hpp"


int main()
{
  using native_mat_t = Eigen::Matrix<double,10,10>;
  using mat_t = core::matrix<native_mat_t>;
  native_mat_t M1;
  mat_t AA;

  // static_assert( std::is_same< mat_t::sc_t, double >::value, "ohhh" );
  // static_assert( mat_t::isEigen == 1, "glddl");
		  
  svd::solver<mat_t, svd::svdKind::EigenJacobi> svdSolve;
  svd::solver<mat_t, svd::svdKind::EigenBDCSVD> svdSolve2;
  
  return 0;
}


// MatrixXf M = MatrixXf::Random(3,2);
// cout << "Here is the matrix m:" << endl << m << endl;
// JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
// cout << "Its singular values are:" << endl << svd.singularValues() << endl;
// cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
// cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
// Vector3f rhs(1, 0, 0);
// cout << "Now consider this rhs vector:" << endl << rhs << endl;
// cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;



// using native_vec_t = Eigen::Vector<double>;
// using native_mat_t = Eigen::Matrix<double>;

// class BG1D{
//   void operator()(const native_vec_t & in, native_vec_t & dudt, double t){
//     // evaluate for a given time 
//   };
//   native_vec_t density_;
// };

// void createSnapshotMatrix(mat_t & myM, BG1D & appObj){
//   // whatever you need to do here
//   auto result = appObj.run();
//   myM
// };

// template <typename app_t, typename mat_t, typename vec_t>
// class galerkinROM{
//   void run(app_t & app, )
//   {
    
//   };
// };

// int main()
// {
  
//   //create the app obj
//   BG1D appObj;

//   // snapshot matrix
//   using mat_t = core::matrix<native_mat_t>;
//   mat_t myM;
//   createSnapshotMatrix(myM, appObj);

//   // svd
//   using svdSolver_t = svd::solver<mat_t>;
//   svdSolver_t svdObj(myM,"Jacobi" or "BDCSVD");
//   svdObj.compute();
//   auto * U = svdObj.leftSingVectors();
//   auto * S = svdObj.singValues();
//   auto * V = svdObj.rightSingVectors();
    
//   return 0;
// };




  // // Eigen::MatrixXd m(2,2);
  // // m(0,0) = 3;
  // // m(1,0) = 2.5;
  // // m(0,1) = -1;
  // // m(1,1) = m(1,0) + m(0,1);
  // // std::cout << m << std::endl;

  // // define the matrix 
  // using wrapped_mat_t = Eigen::Matrix<double,10,10>;
  // wrapped_mat_t AA;
  // using mat_t = core::matrix<wrapped_mat_t>;
  // mat_t myM(AA);

  // // code filling the matrix 
  // // ...

  // // do svd 
  // using svdSolver_t = svd::solver<mat_t>;
  // svdSolver_t svdObj(myM,"Jacobi" or "BDCSVD");
  // svdObj.compute();


  
  
  // using wrap_t = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
  // using mat_t = core::matrix<wrap_t,double>;
  // mat_t A;

  // // using wrap_t2 = vecT<vecT<double>>;
  // // using mat_t2 = core::matrix<wrap_t2,double>;
  // // mat_t2 A2;
  
  // // using wrap_t3 = Eigen::Matrix<double,3,3>;
  // // core::matrix<wrap_t3,double,3,3> DD;
  // // //mat_t3 A3;

  
