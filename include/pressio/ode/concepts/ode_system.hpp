
#ifndef ODE_CONCEPTS_SYSTEM_ALL_HPP_
#define ODE_CONCEPTS_SYSTEM_ALL_HPP_

#include "ode_predicates_for_system.hpp"
#include "ode_has_const_discrete_residual_jacobian_method.hpp"

namespace pressio{ namespace ode{

template<class T, class enable = void>
struct OdeSystem : std::false_type{};

template<class T>
struct OdeSystem<
  T,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct OdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct OdeSystemFusingRhsAndJacobian<
  T,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::rhs_type &>(),
#ifdef PRESSIO_ENABLE_CXX17
		std::declval< std::optional<typename T::jacobian_type*> >()
#else
		std::declval< typename T::jacobian_type* >()
#endif
	       )
	   )
      >::value
   >
  > : std::true_type{};


template<class T, class enable = void>
struct OdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct OdeSystemFusingMassMatrixAndRhs<
  T,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};


template<class T, class enable = void>
struct CompleteOdeSystem : std::false_type{};

template<class T>
struct CompleteOdeSystem<
  T,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>(),
#ifdef PRESSIO_ENABLE_CXX17
		std::declval< std::optional<typename T::jacobian_type*> >()
#else
		std::declval< typename T::jacobian_type* >()
#endif
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, int NumStates,class enable = void>
struct FullyDiscreteSystemWithJacobian : std::false_type{};

template<class T, int NumStates>
struct FullyDiscreteSystemWithJacobian<
  T, NumStates,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && ::pressio::has_discrete_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && std::is_copy_constructible<typename T::discrete_jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_discrete_residual_method_return_result<
      T, typename T::discrete_residual_type>::value
    && ::pressio::ode::has_const_create_discrete_jacobian_method_return_result<
      T, typename T::discrete_jacobian_type>::value
    //
    && ::pressio::ode::has_const_discrete_residual_jacobian_method<
      T, NumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::independent_variable_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      typename T::discrete_jacobian_type
      >::value
    >
  > : std::true_type{};


//
// refine for real-valued case
//
template<class T, class enable = void>
struct RealValuedOdeSystem : std::false_type{};

template<class T>
struct RealValuedOdeSystem<
  T, mpl::enable_if_t<
       OdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedOdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingRhsAndJacobian<
  T,
  mpl::enable_if_t<
    OdeSystemFusingRhsAndJacobian<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  > > : std::true_type{};


template<class T, class enable = void>
struct RealValuedOdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingMassMatrixAndRhs<
  T, mpl::enable_if_t<
    OdeSystemFusingMassMatrixAndRhs<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};


template<class T, class enable = void>
struct RealValuedCompleteOdeSystem : std::false_type{};

template<class T>
struct RealValuedCompleteOdeSystem<
  T, mpl::enable_if_t<
       CompleteOdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

template <class T, int NumStates, class = void>
struct RealValuedFullyDiscreteSystemWithJacobian : std::false_type{};

template <class T, int NumStates>
struct RealValuedFullyDiscreteSystemWithJacobian<
  T, NumStates,
  mpl::enable_if_t<
    FullyDiscreteSystemWithJacobian<T, NumStates>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_residual_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_jacobian_type> >::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};


//
// policy
//
template<class T, class = void>
struct ImplicitResidualJacobianPolicy : std::false_type{};

template<class T>
struct ImplicitResidualJacobianPolicy<
  T,
  ::pressio::mpl::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && all_have_traits_and_same_scalar<
      typename T::state_type,
      typename T::residual_type,
      typename T::jacobian_type>::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type>>::value
    //
    // create methods
    //
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type>::value
    && std::is_same<
      typename T::residual_type,
      decltype(std::declval<T const>().createResidual())
      >::value
    && std::is_same<
      typename T::jacobian_type,
      decltype(std::declval<T const>().createJacobian())
      >::value

    && std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<StepScheme const &>(),
	std::declval<typename T::state_type const &>(),
	std::declval<ImplicitStencilStatesDynamicContainer<typename T::state_type> const & >(),
	std::declval<ImplicitStencilRightHandSideDynamicContainer<typename T::residual_type> & >(),
	std::declval< ::pressio::ode::StepEndAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval<typename T::residual_type &>(),
#ifdef PRESSIO_ENABLE_CXX17
	std::declval< std::optional<typename T::jacobian_type*> >()
#else
	std::declval< typename T::jacobian_type* >()
#endif
	)
       )
      >::value
    >
  > : std::true_type{};

}}
#endif
