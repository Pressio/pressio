
#ifndef PRESSIO_ROM_IMPL_REDUCED_OPERATORS_HELPERS_HPP_
#define PRESSIO_ROM_IMPL_REDUCED_OPERATORS_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class ReducedStateType, class = void>
struct CreateReducedState;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  mpl::enable_if_t< ::pressio::is_vector_eigen<ReducedStateType>::value >
  > {
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType(::pressio::ops::extent(basis, 1));
  }
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  mpl::enable_if_t< ::pressio::is_vector_kokkos<ReducedStateType>::value >
  > {
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType("tmp", ::pressio::ops::extent(basis, 1));
  }
};
#endif
// ------------------------------------------

template<typename T, class = void>
struct galerkin_mass_matrix_type_if_state_is;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename T>
struct galerkin_mass_matrix_type_if_state_is<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif
// ------------------------------------------

template<typename T, class = void>
struct determine_galerkin_jacobian_type_from_state;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename T>
struct determine_galerkin_jacobian_type_from_state<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif
// ------------------------------------------

template<class MassMatrixType, class = void>
struct CreateGalerkinMassMatrix;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class MassMatrixType>
struct CreateGalerkinMassMatrix<
  MassMatrixType,
  mpl::enable_if_t< ::pressio::is_dense_matrix_eigen<MassMatrixType>::value >
  >
{
  template<class BasisType>
  MassMatrixType operator()(const BasisType & basis){
    return MassMatrixType(::pressio::ops::extent(basis, 1),
        ::pressio::ops::extent(basis, 1));
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
  template<class BasisType>
  MassMatrixType operator()(const BasisType & basis){
    return MassMatrixType("tmpMM",
        ::pressio::ops::extent(basis, 1),
        ::pressio::ops::extent(basis, 1));
  }
};
#endif
// ------------------------------------------

template<class RhsType, class = void>
using CreateGalerkinRhs = CreateReducedState<RhsType>;

template<class JacType, class = void>
using CreateGalerkinJacobian = CreateGalerkinMassMatrix<JacType>;


}}} // end pressio::rom::impl
#endif
