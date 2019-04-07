
#ifndef CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{
namespace core{

template<typename derived_type>
class MatrixDenseSharedMemBase
  : private core::details::CrtpBase<
     MatrixDenseSharedMemBase<derived_type>>{
  using this_t = MatrixDenseSharedMemBase<derived_type>;

  static_assert( details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed matrix inheriting from dense sharedMem base!");

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend core::details::CrtpBase<this_t>;

  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;
  MatrixDenseSharedMemBase() = default;
  ~MatrixDenseSharedMemBase() = default;

};//end class

} // end namespace core
}//end namespace rompp
#endif
