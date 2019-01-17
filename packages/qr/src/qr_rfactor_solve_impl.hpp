
#ifndef QR_RFACTOR_SOLVE_HPP_
#define QR_RFACTOR_SOLVE_HPP_

#include "qr_ConfigDefs.hpp"
#include "../../core/src/vector/core_vector_meta.hpp"
#include "../../core/src/matrix/core_matrix_meta.hpp"
#ifdef HAVE_TRILINOS
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include <Teuchos_SerialQRDenseSolver.hpp>
#endif

namespace rompp{ namespace qr{ namespace impl{


template<typename vector_type,
	 typename R_type,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_eigen_vector_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_eigen_dense_matrix_wrapper<R_type>::value
	   > * = nullptr>
void solve(const vector_type & rhs,
	   const R_type & Rmatrix,
	   vector_type & y){
  auto nat_rhs = rhs.data();
  auto nat_y = y.data();
  auto nat_R = Rmatrix.data();

  *nat_y = (*nat_R).template triangularView<Eigen::Upper>().solve(*nat_rhs);
}



#ifdef HAVE_TRILINOS
template<typename vector_type, typename R_type,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_dense_vector_teuchos_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_dense_matrix_teuchos_rcp<R_type>::value
	   > * = nullptr>
void solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  using ord_t	 = typename R_type::element_type::ordinalType;
  using sc_t	 = typename R_type::element_type::scalarType;
  using wrapv_t	 = typename core::details::traits<vector_type>::wrapped_t;
  using solver_t = Teuchos::SerialDenseSolver<ord_t, sc_t>;

  solver_t My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors( Teuchos::rcp( y.data(), false ),
			Teuchos::rcp( const_cast<wrapv_t*>(rhs.data()), false ) );

  int info = My_Solver.factor();
  if (info != 0)
    std::cout << "SerialDenseSolver::factor() returned : "
	      << info << std::endl;

  info = My_Solver.solve();
  if (info != 0)
    std::cout << "SerialDenseSolver::solve() returned : "
	      << info << std::endl;
}
#endif



#ifdef HAVE_TRILINOS
template<typename vector_type,
	 typename R_type,
	 int n,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_dense_vector_teuchos_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_dense_matrix_teuchos_rcp<R_type>::value and
	   n == core::constants::dynamic
	   > * = nullptr>
void solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  // n is not used here, call the one above
  solve(rhs, Rmatrix, y);
}
#endif



#ifdef HAVE_TRILINOS
template<typename vector_type,
	 typename R_type,
	 int n,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_eigen_vector_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_dense_matrix_teuchos_rcp<R_type>::value and
	   n != core::constants::dynamic and n>=1
	   > * = nullptr>
void solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  //using ord_t = typename R_type::element_type::ordinalType;
  using sc_t  = typename R_type::element_type::scalarType;

  //  auto vecSize = rhs.size();
  using eigMat = Eigen::Matrix<sc_t, n, n>;
  core::Matrix<eigMat> eigR( Rmatrix->values() );
  solve(rhs, eigR, y);
}
#endif



#ifdef HAVE_TRILINOS
template<typename vector_type,
	 typename R_type,
	 int n,
	 ::rompp::core::meta::enable_if_t<
	   ::rompp::core::meta::is_eigen_vector_wrapper<vector_type>::value and
	   ::rompp::core::meta::is_dense_matrix_teuchos_rcp<R_type>::value and
	   n == core::constants::dynamic
	   > * = nullptr>
void solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  // todo: maybe here we should use directly backward substitution but it does
  // not seem to be available directly from teuchos

  using ord_t	  = typename R_type::element_type::ordinalType;
  using sc_t	  = typename R_type::element_type::scalarType;
  using tservec_t = Teuchos::SerialDenseVector<ord_t, sc_t>;

  auto vecSize = rhs.size();
  tservec_t rhsTV(Teuchos::View, const_cast<sc_t*>(rhs.data()->data()), vecSize);
  tservec_t yTV(Teuchos::View, y.data()->data(), vecSize);

  Teuchos::SerialQRDenseSolver<ord_t, sc_t> My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors( Teuchos::rcp(&yTV, false), Teuchos::rcp(&rhsTV, false) );
  int info = My_Solver.factor();
  if (info != 0)
    std::cout << "SerialDenseSolver::factor() returned : "
	      << info << std::endl;

  info = My_Solver.solve();
  if (info != 0)
    std::cout << "SerialDenseSolver::solve() returned : "
	      << info << std::endl;
}
#endif



}}}//end namespace rompp::qr::impl
#endif
