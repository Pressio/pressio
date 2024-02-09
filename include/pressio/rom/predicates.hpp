
#ifndef ROM_PREDICATES_HPP_
#define ROM_PREDICATES_HPP_

#include "./impl/ode_has_const_discrete_residual_jacobian_action_method.hpp"

namespace pressio{ namespace rom{

#define PRESSIO_ROM_IMPL_HAS_CONST_CREATE_RETURN_RESULT(NAMEA, NAMEB)	\
template <class T, class = void> \
struct has_const_create_##NAMEA##_return_result \
  : std::false_type{};\
template <class T>\
struct has_const_create_##NAMEA##_return_result<T,\
  std::enable_if_t<\
    std::is_same<\
      decltype(\
	       std::declval<T const>().create##NAMEB()\
	       ),\
      typename T::NAMEA##_type\
      >::value\
    >> : std::true_type{};\

PRESSIO_ROM_IMPL_HAS_CONST_CREATE_RETURN_RESULT(reduced_state,
						ReducedState)
PRESSIO_ROM_IMPL_HAS_CONST_CREATE_RETURN_RESULT(full_state,
						FullState)
// ---------------------------------------------------------------

template <class T, class = void>
struct has_const_map_from_reduced_state_return_void
  : std::false_type{};

template <class T>
struct has_const_map_from_reduced_state_return_void<
  T,
  std::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().mapFromReducedState
	       (
		std::declval<const typename T::reduced_state_type &>(),
		std::declval<typename T::full_state_type &>()
		)
	       )
      >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class = void>
struct has_const_create_full_state_from_reduced_state
  : std::false_type{};

template <class T>
struct has_const_create_full_state_from_reduced_state<
  T,
  std::enable_if_t<
    std::is_same<
      decltype(
	       std::declval<T const>().createFullStateFromReducedState
	       (
		std::declval<const typename T::reduced_state_type &>()
		)
	       ),
      typename T::full_state_type
    >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class ResidualType, class = void>
struct has_const_create_residual_method_return_result
  : std::false_type{};

template <class T, class ResidualType>
struct has_const_create_residual_method_return_result<
  T, ResidualType,
  std::enable_if_t<
    !std::is_void<ResidualType>::value and
    std::is_same<
      ResidualType,
      decltype(
	       std::declval<T const>().createResidual()
	       )
      >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class StateType, class ResidualType, class = void>
struct has_const_residual_method_accept_state_result_return_void
  : std::false_type{};

template <class T, class StateType, class ResidualType>
struct has_const_residual_method_accept_state_result_return_void<
  T, StateType, ResidualType,
  std::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().residual(
					  std::declval<StateType const&>(),
					  std::declval<ResidualType &>()
					  )
	   )
      >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <
  class T,
  class StateType,
  class OperandType,
  class TimeType,
  class ResultType,
  class = void
  >
struct has_const_apply_jacobian_method_accept_state_operand_time_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class OperandType,
  class TimeType,
  class ResultType
  >
struct has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
  T, StateType, OperandType, TimeType, ResultType,
  ::pressio::mpl::void_t<
    decltype
    (
     std::declval<T const>().applyJacobian
     (
      std::declval<StateType const&>(),
      std::declval<OperandType const&>(),
      std::declval<TimeType const &>(),
      std::declval<ResultType &>()
      )
     )
    >
  >: std::true_type{};
// ---------------------------------------------------------------

template <
  class T,
  class StateType,
  class OperandType,
  class TimeType,
  class ResultType,
  class = void
  >
struct has_const_apply_mass_matrix_method_accept_state_operand_time_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class OperandType,
  class TimeType,
  class ResultType
  >
struct has_const_apply_mass_matrix_method_accept_state_operand_time_result_return_void<
  T, StateType, OperandType, TimeType, ResultType,
  ::pressio::mpl::void_t<
    decltype
    (
     std::declval<T const>().applyMassMatrix
     (
      std::declval<StateType const&>(),
      std::declval<OperandType const&>(),
      std::declval<TimeType const &>(),
      std::declval<ResultType &>()
      )
     )
    >
  >: std::true_type{};
// ---------------------------------------------------------------

template <
  class T,
  class OperandType,
  class ResultType,
  class = void
  >
struct has_const_apply_mass_matrix_method_accept_operand_result_return_void
  : std::false_type{};

template <
  class T,
  class OperandType,
  class ResultType
  >
struct has_const_apply_mass_matrix_method_accept_operand_result_return_void<
  T, OperandType, ResultType,
  ::pressio::mpl::void_t<
    decltype
    (
     std::declval<T const>().applyMassMatrix
     (
      std::declval<OperandType const&>(),
      std::declval<ResultType &>()
      )
     )
    >
  >: std::true_type{};
// ---------------------------------------------------------------

template <class T, class OperandType, class = void>
struct has_const_create_result_of_jacobian_action_on
  : std::false_type{};

template <class T, class OperandType>
struct has_const_create_result_of_jacobian_action_on<
  T, OperandType,
  std::enable_if_t<
    !std::is_void<
      decltype
      (
       std::declval<T const>().createResultOfJacobianActionOn
       (
	std::declval<OperandType const &>()
	)
       )
      >::value
    >
  > : std::true_type{};


template <class T, class OperandType, class = void>
struct has_const_create_apply_mass_matrix_result_method_accept_operand_return_result
  : std::false_type{};

template <class T, class OperandType>
struct has_const_create_apply_mass_matrix_result_method_accept_operand_return_result<
  T, OperandType,
  std::enable_if_t<
    !std::is_void<
      decltype
      (
       std::declval<T const>().createResultOfMassMatrixActionOn
       (
	std::declval<OperandType const &>()
	)
       )
      >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class StateType, class IndVarType, class RhsType, class = void>
struct has_const_rhs_method_accept_state_indvar_result_return_void
  : std::false_type{};

template <class T, class StateType, class IndVarType, class RhsType>
struct has_const_rhs_method_accept_state_indvar_result_return_void<
  T, StateType, IndVarType, RhsType,
  std::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T const>().rhs(
					  std::declval<StateType const&>(),
					  std::declval<IndVarType const &>(),
					  std::declval<RhsType &>()
					  )
	   )
      >::value
    >
  > : std::true_type{};

}} // end pressio::rom
#endif  // ROM_PREDICATES_HPP_
