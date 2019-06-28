
#ifndef ALGEBRA_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_
#define ALGEBRA_MATRIX_BASE_MATRIX_DENSE_DISTRIBUTED_BASE_HPP_

#include "../algebra_matrix_traits.hpp"

namespace rompp{
namespace algebra{

template<typename derived_type>
class MatrixDenseDistributedBase
  : private utils::details::CrtpBase<
  MatrixDenseDistributedBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==0,
  "OOPS: non-distributed matrix inheriting from dense distributed base!");
  static_assert( details::traits<derived_type>::is_dense==1,
  "OOPS: non-dense matrix inheriting from dense distributed base!");

  using traits_t = details::traits<derived_type>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;

public:
  void replaceGlobalValue(GO_t globRow, GO_t globCol, sc_t value){
    this->underlying().replaceGlobalValueImpl(globRow, globCol, value);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<
    MatrixDenseDistributedBase<derived_type>>;

  MatrixDenseDistributedBase() = default;
  ~MatrixDenseDistributedBase() = default;

};//end class

} // end namespace algebra
}//end namespace rompp
#endif
