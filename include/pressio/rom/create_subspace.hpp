
#ifndef ROM_CREATE_SUBSPACE_HPP_
#define ROM_CREATE_SUBSPACE_HPP_

#include "./linear_subspace.hpp"
#include "./impl/linear_trial_column_subspace.hpp"

namespace pressio{ namespace rom{

template<
  class ReducedStateType,
  class BasisMatrixType,
  class FullStateType
>
auto create_trial_column_subspace(BasisMatrixType && basisMatrix,
				  FullStateType && offset,
				  bool isAffine)
{
  using basis_matrix_type = mpl::remove_cvref_t<BasisMatrixType>;
  using full_state_type = mpl::remove_cvref_t<FullStateType>;

  static_assert(ValidReducedState<ReducedStateType>::value,
		"Invalid type for the reduced state");
  static_assert( ::pressio::mpl::all_of_t<
		 std::is_copy_constructible, full_state_type, basis_matrix_type>::value,
		"template arguments must be copy constructible");
  static_assert( !::pressio::mpl::all_of_t<
		 std::is_pointer, full_state_type, basis_matrix_type>::value,
		"template arguments cannot be pointers");
  static_assert( !::pressio::mpl::all_of_t<
		 mpl::is_std_shared_ptr, full_state_type, basis_matrix_type>::value,
		"std::unique_ptr are not valid template arguments");
  static_assert(::pressio::all_have_traits_and_same_scalar<
		ReducedStateType, full_state_type, basis_matrix_type>::value,
		"std::unique_ptr are not valid template arguments");

  using ret_t = impl::TrialColumnSubspace<basis_matrix_type, full_state_type, ReducedStateType>;
  return ret_t(std::forward<BasisMatrixType>(basisMatrix),
	       std::forward<FullStateType>(offset),
	       isAffine);
}

}} // end pressio::rom
#endif  // ROM_CREATE_SUBSPACE_HPP_
