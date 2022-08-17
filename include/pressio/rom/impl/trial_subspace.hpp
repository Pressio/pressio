
#ifndef PRESSIO_ROM_TRIAL_SUBSPACES_IMPL_HPP_
#define PRESSIO_ROM_TRIAL_SUBSPACES_IMPL_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class ReducedStateType, class BasisType, class FullStateType>
class PossiblyAffineTrialSubspace
{

public:
  using reduced_state_type = ReducedStateType;
  using basis_type = BasisType;
  using full_state_type = FullStateType;

  PossiblyAffineTrialSubspace() = delete;

  PossiblyAffineTrialSubspace(const BasisType & phi,
			      const FullStateType & shiftVectorIn,
			      bool isAffine)
    : basis_(::pressio::ops::clone(phi)),
      shiftVector_(::pressio::ops::clone(shiftVectorIn)),
      isAffine_(isAffine){}

  PossiblyAffineTrialSubspace(BasisType && phi,
			      FullStateType && shiftVectorIn,
			      bool isAffine)
    : basis_(std::move(phi)),
      shiftVector_(std::move(shiftVectorIn)),
      isAffine_(isAffine){}

  PossiblyAffineTrialSubspace(const BasisType & phi,
			      FullStateType && shiftVectorIn,
			      bool isAffine)
    : basis_(::pressio::ops::clone(phi)),
      shiftVector_(std::move(shiftVectorIn)),
      isAffine_(isAffine){}

  PossiblyAffineTrialSubspace(BasisType && phi,
			      const FullStateType & shiftVectorIn,
			      bool isAffine)
    : basis_(std::move(phi)),
      shiftVector_(::pressio::ops::clone(shiftVectorIn)),
      isAffine_(isAffine){}

  reduced_state_type createReducedState() const{
    auto result = impl::CreateReducedState<ReducedStateType>()(basis_);
    using sc_t = typename ::pressio::Traits<ReducedStateType>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  full_state_type createFullState() const
  {
    // we need to use clone here because full_state_type might
    // NOT have value semantics so we need to ensure a new object
    // is created every time
    auto result = ::pressio::ops::clone(shiftVector_);
    using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  // here we template on reduced state type because it could be an
  // expression not necessarily the reduced state object
  template<class ReducedStateToMap>
  void mapFromReducedState(const ReducedStateToMap & latState,
                          full_state_type & fullState) const
  {
    // always do y = phi*latState
    mapFromReducedStateWithoutShift(latState, fullState);

    if (isAffine_){
      // update full state to account for shift
      using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
      ::pressio::ops::update(fullState, one, shiftVector_, one);
    }
  }

  template<class ReducedStateToMap>
  full_state_type createFullStateFromReducedState(const ReducedStateToMap & latState) const
  {
    auto fomState = this->createFullState();
    this->mapFromReducedState(latState, fomState);
    return fomState;
  }

  const basis_type & viewBasis() const{
    return basis_;
  }

private:
  template<class ReducedStateToMap>
  void mapFromReducedStateWithoutShift(const ReducedStateToMap & latState,
				       full_state_type & fullState) const
  {

    using basis_sc_t = typename ::pressio::Traits<basis_type>::scalar_type;
    using full_state_sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<basis_sc_t>::one();
    constexpr auto beta  = ::pressio::utils::Constants<full_state_sc_t>::zero();
    ::pressio::ops::product(::pressio::nontranspose(), alpha,
			    basis_, latState, beta, fullState);
  }

private:
  const BasisType basis_;
  const FullStateType shiftVector_;
  bool isAffine_ = {};
};

}}}
#endif
