
#ifndef ROM_OPERATORS_HPP_
#define ROM_OPERATORS_HPP_

#include "rom_ConfigDefs.hpp"
#include "CORE_MATRIX"
#include "experimental/rom_matrix_pseudo_inverse.hpp"

namespace rom{
namespace exp{

// this is the operator: (P phi)^+ P 
template <typename phi_t, typename P_t>
class WeightOperator{
public:
  WeightOperator(phi_t const & phi, P_t const & pop)
    : phi_(&phi), pop_(&pop){}

  ~WeightOperator() = default;

  // template <typename vecT,
  // 	    typename std::enable_if<
  // 	      core::details::traits<vecT>::isVector
  // 	      >::type * = nullptr
  // 	    >
  // auto apply(const vecT & b)
  // {
  //   // // P phi
  //   // auto A1 = core::matrixMatrixProduct(*pop_, *phi_);
  //   // std::cout << "P phi: "
	 //   //    << A1.globalRows() << " "
	 //   //    << A1.globalCols() << std::endl;

  //   // // (P phi)^+ => dense matrix
  //   // auto A2 = pseudoInverse(A1);
  //   // std::cout << "(P phi)^+ :"
	 //   //    << A2.globalRows() << " "
	 //   //    << A2.globalCols() << std::endl;

  //   // // P b => vector
  //   // auto b1 = core::matrixVectorProduct(*pop_, b);
  //   // std::cout << "P b :"
	 //   //    << b1.globalSize() << std::endl;

  //   // // (P phi)^+ P b => vector
  //   // return core::matrixVectorProduct(A2, b1);
  // }

  // template <typename T,
  // 	    typename std::enable_if<
  // 	      core::details::traits<T>::isMatrix
  // 	      >::type * = nullptr
  // 	    >
  // auto apply(const T & B)
  // {
  //   // P phi
  //   auto A1 = core::matrixMatrixProduct(*pop_, *phi_);
  //   std::cout << "P phi: "
	 //      << A1.globalRows() << " "
	 //      << A1.globalCols() << std::endl;

  //   // (P phi)^+ => dense matrix
  //   auto A2 = pseudoInverse(A1);
  //   std::cout << "(P phi)^+ :"
	 //      << A2.globalRows() << " "
	 //      << A2.globalCols() << std::endl;

  //   // P B => matrix
  //   auto B1 = core::matrixMatrixProduct(*pop_, B, false, false);
  //   std::cout << "P B :"
	 //      << B1.globalRows() << " "
	 //      << B1.globalCols() << std::endl;

  //   // (P phi)^+ P B => matrix
  //   return core::matrixVectorProduct(A2, B1);
  // }
  
private:
  phi_t const * phi_;
  P_t const * pop_;
};
  
  
}//end namespace exp
}//end namespace rom
#endif 
