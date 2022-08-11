
#ifndef PRESSIO_ROM_SUBSPACES_HPP_
#define PRESSIO_ROM_SUBSPACES_HPP_

#include "./impl/reduced_operators_helpers.hpp"
#include "./impl/trial_subspace.hpp"

namespace pressio{ namespace rom{

/*
  below we use static asserts to check constraints but this is
  not fully correct because constraints should have an impact on the
  overload resolution read this:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  For now we decide to use static asserts to have readable error messages.
  This is kind of justified here, for now, because we have a single function,
  not an overload set to manipulate via sfinae/concepts.
  Later on we can use C++20 concepts to enfore this.
*/

// template<
//   class ReducedStateType,
//   class FullStateType,
//   class BasisType>
// auto create_trial_subspace(BasisType && basis)
// {
//   static_assert(ValidReducedState<ReducedStateType>::value,
// 		"Invalid type for the reduced state");
//   static_assert(std::is_copy_constructible< mpl::remove_cvref_t<BasisType> >::value,
// 		"Basis type must be copy constructible");
//   static_assert(std::is_copy_constructible< FullStateType >::value,
// 		"Full state type must be copy constructible");

//   // use BasisType as template arg below because it carries the correct
//   // qualification to be used inside the class by instanceOrRefWrapper
//   using ret_t = impl::TrialSubspace<ReducedStateType, BasisType, FullStateType>;
//   return ret_t(std::forward<BasisType>(basis));
// }

template<
  class ReducedStateType,
  class FullStateType,
  class BasisType>
auto create_trial_subspace(BasisType && basis,
			   FullStateType && offset,
			   bool isAffine)
{
  static_assert(ValidReducedState<ReducedStateType>::value,
		"Invalid type for the reduced state");
  static_assert(std::is_copy_constructible< mpl::remove_cvref_t<BasisType> >::value,
		"Basis type must be copy constructible");
  static_assert(std::is_copy_constructible< FullStateType >::value,
		"Full state type must be copy constructible");

  // use BasisType and FullStateType as template args below because they carry
  // the correct qualification to be used inside the class by instanceOrRefWrapper
  using ret_t = impl::PossiblyAffineTrialSubspace<ReducedStateType, BasisType, FullStateType>;
  return ret_t(std::forward<BasisType>(basis), offset, isAffine);
}

}} // end pressio::rom
#endif
