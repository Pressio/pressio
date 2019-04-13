
#ifndef SVD_SOLVER_GENERIC_BASE_HPP_
#define SVD_SOLVER_GENERIC_BASE_HPP_

#include "svd_solver_traits.hpp"

namespace rompp{
namespace svd{

template<typename derived_type>
class SolverBase
  : private core::details::CrtpBase<SolverBase<derived_type>>
{

private:
  using matrix_t = typename svd::details::svd_traits<derived_type>::matrix_t;
  using sc_t = typename svd::details::svd_traits<derived_type>::scalar_t;
  using leftSvec_t = typename svd::details::svd_traits<derived_type>::lsv_t;
  using rightSvec_t = typename svd::details::svd_traits<derived_type>::rsv_t;
  using sval_t = typename svd::details::svd_traits<derived_type>::sval_t;

public:

  template<svdType svd_enum_value,
	   typename std::enable_if<
	     svd_enum_value==svdType::truncated
	     >::type * = nullptr>
  void compute(matrix_t & mat, int t, sc_t tol = 1e-12){
    this->underlying().template computeImpl<svd_enum_value>(mat, t, tol);
  }

  const leftSvec_t & cRefLeftSingularVectors() const {
    return this->underlying().cRefLeftSingularVectorsImpl();
  };

  const rightSvec_t & cRefRightSingularVectors() const {
    return this->underlying().cRefRightSingularVectorsImpl();
  };

  const sval_t & singularValues() const{
    return this->underlying().singularValuesImpl();
  };

private:
  SolverBase() = default;
  ~SolverBase() = default;

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend core::details::CrtpBase<SolverBase<derived_type>>;

};//end class

} // end namespace svd
}//end namespace rompp
#endif
