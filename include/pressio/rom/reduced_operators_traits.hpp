
#ifndef PRESSIO_ROM_GALERKIN_REDUCED_OPERATORS_TRAITS_HPP_
#define PRESSIO_ROM_GALERKIN_REDUCED_OPERATORS_TRAITS_HPP_

namespace pressio{ namespace rom{

/*
  steady galerkin
*/
template<class T, class = void>
struct SteadyGalerkinDefaultOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_residual_type = void;
  using reduced_jacobian_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct SteadyGalerkinDefaultOperatorsTraits<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_residual_type = T;

  // figure out what is the reduced jacobian type
  // if the reduced state is Eigen vector,
  // it makes sense to use an Eigen dense matrix to store
  // the Galerkin jacobian since all reduced operators are dense
  using reduced_jacobian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

/*
  steady LSPG
*/
template<class T, class = void>
struct SteadyLspgDefaultOperatorsTraits
{
  using reduced_state_type = void;
  using hessian_type    = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct SteadyLspgDefaultOperatorsTraits<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type = T;
  using gradient_type = T;
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif


}}
#endif
