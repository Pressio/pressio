
#ifndef PRESSIO_ROM_SUBSPACES_HPP_
#define PRESSIO_ROM_SUBSPACES_HPP_

#include "./impl/reduced_operators_helpers.hpp"
#include "./impl/trial_subspace.hpp"

namespace pressio{ namespace rom{

/*
  below we kind of abuse notation a bit, since we use static asserts
  also to check constraints. This is not fully correct because
  constraints should have an impact on the overload resolution read this:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  For now we do this to have readable error messages.
  Also, this is kind of justified here, for now, because we have
  a single function, not an overload set to manipulate via sfinae/concepts.
  Later on we can use C++20 concepts to enfore this.
*/

template<
  class ReducedStateType,
  class BasisType,
  class FullStateType
>
auto create_trial_subspace(BasisType && basis,
			   FullStateType && offset,
			   bool isAffine)
{
  // we need to use decay here to get the type after removing
  // reference and cv qualification
  using basis_type = std::decay_t<BasisType>;
  using full_state_type = std::decay_t<FullStateType>;

  // constraints
  static_assert(ValidReducedState<ReducedStateType>::value,
		"Invalid type for the reduced state");
  static_assert(std::is_copy_constructible< basis_type >::value,
		"Basis type must be copy constructible");
  static_assert(std::is_copy_constructible< full_state_type >::value,
		"Full state type must be copy constructible");

  // mandates
  static_assert(std::is_same< typename pressio::Traits<basis_type>::scalar_type,
		typename pressio::Traits<full_state_type>::scalar_type >::value,
		"Mismatching scalar_type");

  using ret_t = impl::PossiblyAffineTrialSubspace<ReducedStateType, basis_type, full_state_type>;
  return ret_t(std::forward<BasisType>(basis),
	       std::forward<FullStateType>(offset),
	       isAffine);
}

}} // end pressio::rom
#endif
