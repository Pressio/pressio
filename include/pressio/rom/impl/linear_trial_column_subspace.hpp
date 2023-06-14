
#ifndef ROM_IMPL_LINEAR_TRIAL_COLUMN_SUBSPACE_HPP_
#define ROM_IMPL_LINEAR_TRIAL_COLUMN_SUBSPACE_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class ReducedStateType, class = void>
struct CreateReducedState;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  mpl::enable_if_t< ::pressio::is_vector_eigen<ReducedStateType>::value >
  >
{
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType(::pressio::ops::extent(basis, 1));
  }
};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<class ReducedStateType>
struct CreateReducedState<
  ReducedStateType,
  mpl::enable_if_t< ::pressio::is_vector_kokkos<ReducedStateType>::value >
  >
{
  template<class BasisType>
  ReducedStateType operator()(const BasisType & basis){
    return ReducedStateType("tmp", ::pressio::ops::extent(basis, 1));
  }
};
#endif


template <class BasisMatrixType, class FullStateType, class ReducedStateType>
class TrialColumnSubspace
{
public:
  using reduced_state_type = ReducedStateType;
  using basis_matrix_type  = std::remove_cv_t<BasisMatrixType>;
  using full_state_type    = std::remove_cv_t<FullStateType>;

private:
  using linear_subspace_t = LinearSubspace<basis_matrix_type>;
  const linear_subspace_t linSpace_;
  const full_state_type translation_;
  bool isAffine_;
  basis_matrix_type * dummy_ = nullptr;

public:
  TrialColumnSubspace(const basis_matrix_type & basisMatrix,
		      const full_state_type & translation,
		      bool isAffine)
    : linSpace_(basisMatrix,
		linear_subspace_t::SpanningSet::Columns),
      translation_(translation),
      isAffine_(isAffine){}

  TrialColumnSubspace(basis_matrix_type && basisMatrix,
		      full_state_type && translation,
		      bool isAffine)
    : linSpace_(std::move(basisMatrix),
		linear_subspace_t::SpanningSet::Columns),
      translation_(std::move(translation)),
      isAffine_(isAffine){}

  TrialColumnSubspace(const basis_matrix_type & basisMatrix,
		      full_state_type && translation,
		      bool isAffine)
    : linSpace_(basisMatrix,
		linear_subspace_t::SpanningSet::Columns),
      translation_(std::move(translation)),
      isAffine_(isAffine){}

  TrialColumnSubspace(basis_matrix_type && basisMatrix,
		      const full_state_type & translation,
		      bool isAffine)
    : linSpace_(std::move(basisMatrix),
		linear_subspace_t::SpanningSet::Columns),
      translation_(translation),
      isAffine_(isAffine){}

  TrialColumnSubspace(const TrialColumnSubspace & other)
    : linSpace_(other.linSpace_),
      translation_(::pressio::ops::clone(other.translation_)),
      isAffine_(other.isAffine_){}

  /* For this class to be really immutable, we should not have a
     move assign operator and copy assign.
     Since we have a const member, the compiler defines all those as deleted.
     And since we have a copy constructor, the move constructor does not
     particupare in OR so the copy constructor is always called.
  */

  ~TrialColumnSubspace() = default;

  //
  // methods
  //
  reduced_state_type createReducedState() const{
    const auto & basis = linSpace_.basis();
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
    auto result = ::pressio::ops::clone(translation_);
    using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    ::pressio::ops::fill(result, sc_t(0));
    return result;
  }

  void mapFromReducedState(const reduced_state_type & latState,
			   full_state_type & fullState) const
  {
    // always do y = phi*latState
    mapFromReducedStateWithoutTranslation(latState, fullState);

    if (isAffine_){
      // update full state to account for translation
      using sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
      ::pressio::ops::update(fullState, one, translation_, one);
    }
  }

  full_state_type createFullStateFromReducedState(const reduced_state_type & latState) const
  {
    auto fomState = this->createFullState();
    this->mapFromReducedState(latState, fomState);
    return fomState;
  }

  bool isColumnSpace() const{ return true; }
  bool isRowSpace() const{ return false; }
  const full_state_type & translationVector() const{ return translation_; }
  std::size_t dimension() const{ return linSpace_.dimension(); }

  const basis_matrix_type & basisOfTranslatedSpace() const{
    return linSpace_.basis();
  }

  const basis_matrix_type & basis() const{
    if (isAffine_){
      return linSpace_.basis();
    } else{
      return *dummy_;
    }
  }

private:
  void mapFromReducedStateWithoutTranslation(const reduced_state_type & latState,
					     full_state_type & fullState) const
  {

    const auto & basis = linSpace_.basis();
    using basis_sc_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    using full_state_sc_t = typename ::pressio::Traits<full_state_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<basis_sc_t>::one();
    constexpr auto beta  = ::pressio::utils::Constants<full_state_sc_t>::zero();
    ::pressio::ops::product(::pressio::nontranspose(), alpha,
			    basis, latState, beta, fullState);
  }
};

}}}
#endif  // ROM_IMPL_LINEAR_TRIAL_COLUMN_SUBSPACE_HPP_
