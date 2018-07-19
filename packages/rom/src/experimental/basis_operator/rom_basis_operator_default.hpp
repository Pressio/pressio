
#ifndef ROM_BASIS_OPERATOR_DEFAULT_HPP_
#define ROM_BASIS_OPERATOR_DEFAULT_HPP_

#include "rom_ConfigDefs.hpp"
#include "./rom_basis_operator_base.hpp"

namespace rom{
namespace experimental{

template<typename matrix_type,
	  typename
	  std::enable_if<
	    core::meta::is_coreMatrixWrapper<matrix_type>::value && 
	    core::details::traits<matrix_type>::isEigen
	    >::type * = nullptr
	  >
class basisOperatorDefault
  : public basisOperatorBase<basisOperatorDefault<matrix_type>>
{

public:
  basisOperatorDefault(const matrix_type & phi)
    : phi_(&phi), phiT_(phi.cols(), phi.rows()){
    *phiT_.data() = phi_->data()->transpose();
  }
  ~basisOperatorDefault() = default;

private:

  template <typename T>
  void projectImpl(const T & yin, T & yout) const
  {
    // assert( yin.size() == phiT_.cols() );
    // if (yout.size() != phiT_.rows() )
    //   yout.resize(phiT_.rows());

    // project impl
    auto tmp = (*phiT_.data()) * (*yin.data());
    (*yout.data()) = tmp.sparseView();
  }

  template <typename T>
  void leftMultiplyImpl(const T & yin, T & yout) 
  {
    // assert( yin.size() == phi_->rows() );
    // if (yout.size() != phi_->rows() )
    //   yout.resize(phi_->rows());
	
    // multiply impl
    (*yout.data()) = (*phi_->data()) * (*yin.data());
  }

  template <typename T, typename T1>
  void rightMultiplyImpl(const T & yin, T1 & yout) 
  {
    // multiply impl
    auto tmp = (*yin.data()) * (*phi_->data());
    (*yout.data()) = tmp.sparseView();
  }

  
private:
  matrix_type const * phi_;
  matrix_type phiT_;
  
private:
  friend basisOperatorBase<basisOperatorDefault<matrix_type>>;

};//end class

}//end namespace exp
}//end namespace rom

#endif
