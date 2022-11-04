
#ifndef ROM_REDUCED_OPERATORS_TRAITS_HPP_
#define ROM_REDUCED_OPERATORS_TRAITS_HPP_

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

namespace impl{
template<class SubspaceType>
using steady_galerkin_default_reduced_state_t =
  typename SteadyGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_state_type;

template<class SubspaceType>
using steady_galerkin_default_reduced_residual_t =
  typename SteadyGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_residual_type;

template<class SubspaceType>
using steady_galerkin_default_reduced_jacobian_t =
  typename SteadyGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_jacobian_type;
} //end namespace pressio:rom::impl


/*
  unsteady explicit galerkin
*/
template<class T, class = void>
struct ExplicitGalerkinDefaultOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_right_hand_side_type = void;
  using reduced_mass_matrix_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ExplicitGalerkinDefaultOperatorsTraits<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_right_hand_side_type = T;
  using reduced_mass_matrix_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

namespace impl{
template<class SubspaceType>
using explicit_galerkin_default_reduced_state_t =
  typename ExplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_state_type;

template<class SubspaceType>
using explicit_galerkin_default_reduced_right_hand_side_t =
  typename ExplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_right_hand_side_type;

template<class SubspaceType>
using explicit_galerkin_default_reduced_mass_matrix_t =
  typename ExplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_mass_matrix_type;
} //end namespace pressio:rom::impl


/*
  unsteady implicit galerkin
*/
template<class T, class = void>
struct ImplicitGalerkinDefaultOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_residual_type = void;
  using reduced_jacobian_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ImplicitGalerkinDefaultOperatorsTraits<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_residual_type = T;
  using reduced_jacobian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

namespace impl{
template<class SubspaceType>
using implicit_galerkin_default_reduced_state_t =
  typename ImplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_state_type;

template<class SubspaceType>
using implicit_galerkin_default_reduced_residual_t =
  typename ImplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_residual_type;

template<class SubspaceType>
using implicit_galerkin_default_reduced_jacobian_t =
  typename ImplicitGalerkinDefaultOperatorsTraits<
    typename SubspaceType::reduced_state_type>::reduced_jacobian_type;
} //end namespace pressio:rom::impl


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

/*
  unsteady LSPG
*/
template<class T, class = void>
struct UnsteadyLspgDefaultOperatorsTraits
{
  using reduced_state_type = void;
  using hessian_type    = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct UnsteadyLspgDefaultOperatorsTraits<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type = T;
  using gradient_type = T;
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

}}
#endif  // ROM_REDUCED_OPERATORS_TRAITS_HPP_
