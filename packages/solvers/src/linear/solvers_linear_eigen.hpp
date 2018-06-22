
#ifndef SOLVERS_LINEAR_EIGEN_HPP_
#define SOLVERS_LINEAR_EIGEN_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_linear_base.hpp"
#include "matrix/core_matrix_traits.hpp"
#include "vector/core_vector_traits.hpp"

#include<Eigen/SparseLU>
#include<Eigen/SparseCholesky>


namespace solvers {
namespace experimental{
  

//***********************************************
// for DENSE matrix
//***********************************************
template<typename matrix_type, 
	 typename rhs_type, 
	 typename result_type,
	 int algo_type
	 >
class linearSolver<matrix_type, rhs_type, result_type, algo_type,
		   typename 
		   std::enable_if<
		     core::details::traits<matrix_type>::isMatrix==1 &&
		     core::details::traits<matrix_type>::isEigen==1 &&
		     core::details::traits<matrix_type>::isDense==1 &&
		     core::details::traits<rhs_type>::isVector==1 &&
		     core::details::traits<rhs_type>::isEigen==1 &&
		     core::details::traits<result_type>::isVector==1 &&
		     core::details::traits<result_type>::isEigen==1
		     >::type
		   >
  : public linearBase<linearSolver<matrix_type,rhs_type,result_type,algo_type>,
		      matrix_type, rhs_type, result_type>
{
public:
  linearSolver() = default;
  ~linearSolver() = default;

private:
  void solveImpl(const matrix_type & A,
		 const rhs_type & b,
		 result_type & x)const{
    *x.data() = A.data()->colPivHouseholderQr().solve(*b.data());
  }
  friend linearBase<linearSolver<matrix_type,rhs_type,result_type,algo_type>,
		    matrix_type, rhs_type, result_type>;
};

  
//***********************************************
// for SPARSE matrix
//***********************************************
template<typename matrix_type, 
	 typename rhs_type, 
	 typename result_type,
	 int algo_type>
class linearSolver<matrix_type, rhs_type, result_type, algo_type,
		   typename 
		   std::enable_if<
		     core::details::traits<matrix_type>::isMatrix==1 &&
		     core::details::traits<matrix_type>::isEigen==1 &&
		     core::details::traits<matrix_type>::isSparse==1 &&
		     core::details::traits<rhs_type>::isVector==1 &&
		     core::details::traits<rhs_type>::isEigen==1 &&
		     core::details::traits<result_type>::isVector==1 &&
		     core::details::traits<result_type>::isEigen==1
		     >::type
		   >		   
  : public linearBase<linearSolver<matrix_type,rhs_type,
				   result_type, algo_type>,
		      matrix_type, rhs_type, result_type>
{
  using scalar_t = typename core::details::traits<matrix_type>::scalar_t;
  using ordinal_t = typename core::details::traits<matrix_type>::ordinal_t;

  using eigMat_t = typename std::conditional<
    core::details::traits<matrix_type>::isRowMajor==1,
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, ordinal_t>,
    Eigen::SparseMatrix<scalar_t, Eigen::ColMajor, ordinal_t>
    >::type;
public:
  linearSolver() = default;
  ~linearSolver() = default;

private:
  void solveImpl(const matrix_type & A,
		 const rhs_type & b,
		 result_type & x) const
  {
    // //SparseLU<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index>> solver;
    Eigen::ConjugateGradient<eigMat_t, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(*A.data());
 
    //Eigen::BiCGSTAB<eigMat_t> solver;
    //Eigen::SimplicialLDLT<eigMat_t> solver;
    //solver.compute(*A.data());
    *x.data() = solver.solve(*b.data());
  }
  friend linearBase<linearSolver<matrix_type,rhs_type,
				 result_type, algo_type>,
		    matrix_type, rhs_type, result_type>;
};
//**************************************************************


} // end namespace experimental
} // end namespace solvers
#endif
