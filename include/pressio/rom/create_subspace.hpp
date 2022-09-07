
#ifndef PRESSIO_ROM_CREATE_SUBSPACE_HPP_
#define PRESSIO_ROM_CREATE_SUBSPACE_HPP_

#include "./impl/reduced_operators_helpers.hpp"
#include "./linear_subspace.hpp"
#include "./linear_trial_column_subspace.hpp"

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
auto create_trial_column_subspace(BasisType && basis,
				  FullStateType && offset,
				  bool isAffine)
{
  using basis_type         = mpl::remove_cvref_t<BasisType>;
  using full_state_type    = mpl::remove_cvref_t<FullStateType>;
  using ret_t = TrialColumnSubspace<ReducedStateType, basis_type, full_state_type>;
  return ret_t(std::forward<BasisType>(basis),
	       std::forward<FullStateType>(offset),
	       isAffine);
}

}} // end pressio::rom
#endif
