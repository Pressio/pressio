
#ifndef PRESSIO_ROM_SUBSPACES_HPP_
#define PRESSIO_ROM_SUBSPACES_HPP_

#include "./impl/reduced_operators_helpers.hpp"
#include "./impl/trial_subspace.hpp"

namespace pressio{ namespace rom{

template<
  class ReducedStateType, class FullStateType, class BasisType,
  mpl::enable_if_t<
       impl::valid_reduced_state_type<ReducedStateType>::value
    && std::is_copy_constructible< mpl::remove_cvref_t<BasisType> >::value
    && std::is_copy_constructible< FullStateType >::value,
    int > = 0
  >
auto create_trial_subspace(BasisType && basis)
{
  // here we need to use BasisType because it carries the correct
  // qualification to be used inside the class by instanceOrRefWrapper

  using ret_t = impl::TrialSubspace<ReducedStateType, BasisType, FullStateType>;
  return ret_t(std::forward<BasisType>(basis));
}

template<
  class ReducedStateType, class FullStateType, class BasisType,
  mpl::enable_if_t<
       impl::valid_reduced_state_type<ReducedStateType>::value
    && std::is_copy_constructible< mpl::remove_cvref_t<BasisType> >::value
    && std::is_copy_constructible< FullStateType >::value,
    int > = 0
  >
auto create_affine_trial_subspace(BasisType && basis, const FullStateType & shift)
{
  // here we need to use BasisType because it carries the correct
  // qualification to be used inside the class by instanceOrRefWrapper

  using ret_t = impl::AffineTrialSubspace<ReducedStateType, BasisType, FullStateType>;
  return ret_t(std::forward<BasisType>(basis), shift);
}

}} // end pressio::rom
#endif
