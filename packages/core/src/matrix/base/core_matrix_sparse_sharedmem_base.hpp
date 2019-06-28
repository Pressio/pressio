
#ifndef CORE_MATRIX_BASE_MATRIX_SPARSE_SHAREDMEM_BASE_HPP_
#define CORE_MATRIX_BASE_MATRIX_SPARSE_SHAREDMEM_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{
namespace core{

template<typename derived_type>
class MatrixSparseSharedMemBase
  : private utils::details::CrtpBase<
  MatrixSparseSharedMemBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed matrix inheriting from sparse sharedMem base!");

  using traits_t = details::traits<derived_type>;
  using ord_t = typename traits_t::ordinal_t;
  using sc_t = typename traits_t::scalar_t;

public:

  bool isCompressed() const{
    return this->underlying().isCompressedImpl();
  }

  void compress(){
    this->underlying().compressImpl();
  }

  //-----------------------------------------------------------
  // note this insert one by one might not be the best method
  // for efficiency. But it provides a simple nice way to store.
  // NOTE: targetLocation can be either a row index or a columnm
  // depending on whether the matrix is stored row-wise of columnwise.
  //-----------------------------------------------------------
  void insertValues(ord_t targetLocation,
		    ord_t numEntries,
		    const sc_t * values,
		    const ord_t * indices){
    this->underlying().insertValuesImpl(targetLocation,
					numEntries,
					values,
					indices);
  }

  // NOTE: we return by copy. We do not enable reference []
  // because it makes little sense for a sparse matrix
  sc_t operator() (ord_t row, ord_t col) const;

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<
    MatrixSparseSharedMemBase<derived_type>>;

  MatrixSparseSharedMemBase() = default;
  ~MatrixSparseSharedMemBase() = default;

};//end class

} // end namespace core
}//end namespace rompp
#endif
