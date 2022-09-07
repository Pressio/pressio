
#ifndef PRESSIO_ROM_TRIAL_COLUMN_SUBSPACE_HPP_
#define PRESSIO_ROM_TRIAL_COLUMN_SUBSPACE_HPP_

namespace pressio{ namespace rom{

template <class ReducedStateType, class BasisType, class FullStateType>
class TrialColumnSubspace
  : public LinearAffineSubspace<BasisType, FullStateType>
{
  // constraints
  static_assert(ValidReducedState<ReducedStateType>::value,
		"Invalid type for the reduced state");

  using base_t = LinearAffineSubspace<BasisType, FullStateType>;

public:
  using typename base_t::basis_type;
  using typename base_t::offset_type;
  using reduced_state_type = ReducedStateType;
  using full_state_type = std::remove_cv_t<FullStateType>;

  TrialColumnSubspace() = delete;

  TrialColumnSubspace(const basis_type & basis,
		      const full_state_type & offset,
		      bool isAffine)
    : base_t(basis, offset, isAffine){}

  TrialColumnSubspace(basis_type && basis,
		      full_state_type && offset,
		      bool isAffine)
    : base_t(std::move(basis), std::move(offset), isAffine){}

  TrialColumnSubspace(const basis_type & basis,
		      full_state_type && offset,
		      bool isAffine)
    : base_t(basis, std::move(offset), isAffine){}

  TrialColumnSubspace(basis_type && basis,
		      const full_state_type & offset,
		      bool isAffine)
    : base_t(std::move(basis), offset, isAffine){}

  reduced_state_type createReducedState() const{
    const auto & basis = base_t::viewBasis();
    auto result = impl::CreateReducedState<ReducedStateType>()(basis);
    using sc_t = typename ::pressio::Traits<ReducedStateType>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  full_state_type createFullState() const
  {
    // we need to use clone here because full_state_type might
    // NOT have value semantics so we need to ensure a new object
    // is created every time
    const auto & offset = base_t::viewOffset();
    auto result = ::pressio::ops::clone(offset);
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
    mapFromReducedStateWithoutOffset(latState, fullState);

    if (base_t::isAffine()){
      // update full state to account for offset
      using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
      const auto & offset = base_t::viewOffset();
      ::pressio::ops::update(fullState, one, offset, one);
    }
  }

  template<class ReducedStateToMap>
  full_state_type createFullStateFromReducedState(const ReducedStateToMap & latState) const
  {
    auto fomState = this->createFullState();
    this->mapFromReducedState(latState, fomState);
    return fomState;
  }

private:
  template<class ReducedStateToMap>
  void mapFromReducedStateWithoutOffset(const ReducedStateToMap & latState,
					full_state_type & fullState) const
  {

    const auto & basis = base_t::viewBasis();
    using basis_sc_t = typename ::pressio::Traits<basis_type>::scalar_type;
    using full_state_sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<basis_sc_t>::one();
    constexpr auto beta  = ::pressio::utils::Constants<full_state_sc_t>::zero();
    ::pressio::ops::product(::pressio::nontranspose(), alpha,
			    basis, latState, beta, fullState);
  }
};

}}
#endif
