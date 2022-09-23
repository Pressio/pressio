
#ifndef PRESSIO_ROM_IMPL_REDUCED_OPERATORS_HELPERS_HPP_
#define PRESSIO_ROM_IMPL_REDUCED_OPERATORS_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

// --------------------------------------------------------------
// CreateGalerkinRhs
// --------------------------------------------------------------
template<class T, class = void> struct CreateGalerkinRhs;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct CreateGalerkinRhs<
  T, mpl::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  >
{
  T operator()(std::size_t ext){ return T(ext); }
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<class T>
struct CreateGalerkinRhs<
  T, mpl::enable_if_t< ::pressio::is_vector_kokkos<T>::value >
  >
{
  ReducedStateType operator()(std::size_t ext){
    return ReducedStateType("tmp", ext);
  }
};
#endif

// --------------------------------------------------------------
// CreateGalerkinMassMatrix
// --------------------------------------------------------------
template<class MassMatrixType, class = void>
struct CreateGalerkinMassMatrix;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class MassMatrixType>
struct CreateGalerkinMassMatrix<
  MassMatrixType,
  mpl::enable_if_t< ::pressio::is_dense_matrix_eigen<MassMatrixType>::value >
  >
{
  MassMatrixType operator()(std::size_t ext){
    return MassMatrixType(ext, ext);
  }
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<class MassMatrixType>
struct CreateGalerkinMassMatrix<
  MassMatrixType,
  mpl::enable_if_t< ::pressio::is_dense_matrix_kokkos<MassMatrixType>::value >
  >
{
  MassMatrixType operator()(std::size_t ext){
    return MassMatrixType("tmpMM", ext, ext);
  }
};
#endif
// ------------------------------------------

// this is an alias becuase it can done in the same way
template<class JacType, class = void>
using CreateGalerkinJacobian = CreateGalerkinMassMatrix<JacType>;

}}} // end pressio::rom::impl
#endif
