
#ifndef PRESSIO_ROM_PREDICATES_HPP_
#define PRESSIO_ROM_PREDICATES_HPP_

namespace pressio{ namespace rom{

#define PRESSIO_ROM_IMPL_HAS_CONST_CREATE_RETURN_RESULT(NAMEA, NAMEB)	\
template <class T, class = void> \
struct has_const_create_##NAMEA##_return_result \
  : std::false_type{};\
template <class T>\
struct has_const_create_##NAMEA##_return_result<T,\
  mpl::enable_if_t<\
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
  mpl::enable_if_t<
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
  mpl::enable_if_t<
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

template <class T, class = void>
struct has_const_view_basis
  : std::false_type{};

template <class T>
struct has_const_view_basis<
  T,
  mpl::enable_if_t<
    std::is_same<
      decltype(
	       std::declval<T const>().viewBasis()
	       ),
      const typename T::basis_type &
    >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class = void>
struct has_const_view_affine_offset : std::false_type{};

template <class T>
struct has_const_view_affine_offset<
  T,
  mpl::enable_if_t<
    std::is_same<
      decltype(
	       std::declval<T const>().viewAffineOffset()
	       ),
      const typename T::full_state_type &
    >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class = void>
struct has_update_and_view_manifold_tangent_accept_state_return_ref
  : std::false_type{};

template <class T>
struct has_update_and_view_manifold_tangent_accept_state_return_ref<
  T,
  mpl::enable_if_t<
    std::is_same<
      decltype(
	       std::declval<T const>().updateAndViewManifoldTangent
	       (
		std::declval<const typename T::latent_state_type &>()
		)
	       ),
      const typename T::manifold_tangent_type &
    >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------

template <class T, class = void>
struct has_view_manifold_tangent_return_ref
  : std::false_type{};

template <class T>
struct has_view_manifold_tangent_return_ref<
  T,
  mpl::enable_if_t<
    std::is_same<
      decltype(
	       std::declval<T const>().viewManifoldTangent()
	       ),
      const typename T::manifold_tangent_type &
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
  ::pressio::mpl::enable_if_t<
    !std::is_void<ResidualType>::value and
    mpl::is_same<
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
  ::pressio::mpl::enable_if_t<
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
  class ResultType,
  class = void
  >
struct has_const_apply_jacobian_method_accept_state_operand_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class OperandType,
  class ResultType
  >
struct has_const_apply_jacobian_method_accept_state_operand_result_return_void<
  T, StateType, OperandType, ResultType,
  ::pressio::mpl::void_t<
    decltype
    (
     std::declval<T const>().applyJacobian
     (
      std::declval<StateType const&>(),
      std::declval<OperandType const&>(),
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
struct has_const_create_apply_jacobian_result_method_accept_operand_return_result
  : std::false_type{};

template <class T, class OperandType>
struct has_const_create_apply_jacobian_result_method_accept_operand_return_result<
  T, OperandType,
  mpl::enable_if_t<
    !std::is_void<
      decltype
      (
       std::declval<T const>().createApplyJacobianResult
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
  mpl::enable_if_t<
    !std::is_void<
      decltype
      (
       std::declval<T const>().createApplyMassMatrixResult
       (
	std::declval<OperandType const &>()
	)
       )
      >::value
    >
  > : std::true_type{};
// ---------------------------------------------------------------


template <
  class T,
  int n,
  class StepType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class ManifoldJacobian,
  class JacobianType,
  class = void
  >
struct has_const_discrete_residual_jacobian_action_method : std::false_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 1, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteTimeResidualAndJacobianAction
       (
	std::declval<StepType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<ResidualType &>(),
	std::declval<ManifoldJacobian const &>(),
	std::declval<bool>(),
	std::declval<JacobianType &>(),
	std::declval<StateType const&>()
	)
       )
      >::value
    >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 2, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval<bool>(),
       std::declval<JacobianType &>(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
      >::value
  >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 3, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval<bool>(),
       std::declval<JacobianType &>(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
  >::value
  >
  > : std::true_type{};


}} // end pressio::rom
#endif
