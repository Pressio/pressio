
#ifndef PRESSIO_ROM_TRIAL_SUBSPACES_IMPL_HPP_
#define PRESSIO_ROM_TRIAL_SUBSPACES_IMPL_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class ReducedStateType, class BasisType, class FullStateType>
class TrialSubspace
{
  ::pressio::utils::InstanceOrReferenceWrapper<BasisType> basis_;
  
public:
  using reduced_state_type = ReducedStateType;
  using basis_type = mpl::remove_cvref_t<BasisType>;
  using full_state_type = FullStateType;

  explicit TrialSubspace(BasisType && phi)
    : basis_(std::forward<BasisType>(phi)){}

  reduced_state_type createReducedState() const{
    auto result = impl::CreateReducedState<ReducedStateType>()(basis_.get());
    using sc_t = typename ::pressio::Traits<ReducedStateType>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  full_state_type createFullState() const
  {
    using sc_t = typename ::pressio::Traits<reduced_state_type>::scalar_type;
    constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
    auto tmp = impl::CreateReducedState<ReducedStateType>()(basis_.get());
    auto result = ::pressio::ops::product<FullStateType>(::pressio::nontranspose(),
							 one, basis_.get(), tmp);
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  // here we template on reduced state type because it could be an
  // expression not necessarily the reduced state object
  template<class ReducedStateToMap>
  void mapFromReducedState(const ReducedStateToMap & latState,
                          full_state_type & fullState) const
  {
    // NOT affine, do y = phi*latState

    using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), one,
			    basis_.get(), latState, zero, fullState);
  }

  template<class ReducedStateToMap>
  full_state_type createFullStateFromReducedState(const ReducedStateToMap & latState) const
  {
    auto fomState = createFullState();
    mapFromReducedState(latState, fomState);
    return fomState;
  }

  const basis_type & viewBasis() const{
    return basis_.get();
  }
};
// ------------------------------------------


template <class ReducedStateType, class BasisType, class FullStateType>
class AffineTrialSubspace
  : public TrialSubspace<ReducedStateType, BasisType, mpl::remove_cvref_t<FullStateType>>
{
  using base_t = TrialSubspace<ReducedStateType, BasisType, mpl::remove_cvref_t<FullStateType>>;
  ::pressio::utils::InstanceOrReferenceWrapper<FullStateType> shiftVector_;

public:
  using typename base_t::reduced_state_type;
  using typename base_t::basis_type;
  using full_state_type = mpl::remove_cvref_t<FullStateType>;

  AffineTrialSubspace(BasisType && phi,
		      FullStateType && shiftVectorIn)
    : base_t(std::forward<BasisType>(phi)),
      shiftVector_(std::forward<FullStateType>(shiftVectorIn)){}

  full_state_type createFullState() const
  {
    // we need to use clone here because full_state_type might
    // NOT have value semantics so we need to ensure a new object
    // is created every time
    auto result = ::pressio::ops::clone(shiftVector_.get());
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

    // do y = phi*latState + shift
    using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
    base_t::mapFromReducedState(latState, fullState);
    ::pressio::ops::update(fullState, one, shiftVector_.get(), one);
  }

  template<class ReducedStateToMap>
  full_state_type createFullStateFromReducedState(const ReducedStateToMap & latState) const
  {
    auto fomState = this->createFullState();
    this->mapFromReducedState(latState, fomState);
    return fomState;
  }
};

}}}
#endif
