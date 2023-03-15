
#ifndef SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_

#include <typeinfo>
#include <optional>

namespace pressio{
namespace nonlinearsolvers{

template<class T, class = void>
struct normal_eqs_default_types{
  using hessian_type  = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct normal_eqs_default_types<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using gradient_type = T;

  static hessian_type createHessian(const T & v){
    const auto ext = ::pressio::ops::extent(v, 0);
    return hessian_type(ext, ext);
  }
};
#endif

template<class T> using normal_eqs_default_hessian_t =
  typename normal_eqs_default_types<T>::hessian_type;
template<class T> using normal_eqs_default_gradient_t =
  typename normal_eqs_default_types<T>::gradient_type;

// ====================================================================
// ====================================================================
// ====================================================================

template<class T, class = void>
struct valid_state_for_least_squares : std::false_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct valid_state_for_least_squares<
  T, mpl::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  > : std::true_type{};
#endif


// ====================================================================
// ====================================================================
// ====================================================================

namespace impl{

template<class, class> struct TagBasedStaticRegistry;

template<class ...Tags, class ...DataTypes>
struct TagBasedStaticRegistry< std::tuple<Tags...>, std::tuple<DataTypes...> >
{
  std::tuple<DataTypes...> d_;

  template<class ...CArgs>
  TagBasedStaticRegistry(CArgs && ... cargs) : d_(std::forward<CArgs>(cargs)...){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same, Tags...>::value)
      < mpl::size<Tags...>::value;
  }

  template<class TagToFind, class T>
  void set(T && o){
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    std::get<i>(d_) = std::forward<T>(o);
  }

  template<class TagToFind>
  auto & get(){
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    return std::get<i>(d_);
  }

  template<class TagToFind>
  const auto & get() const {
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    return std::get<i>(d_);
  }
};

template<class, class, class> struct TagBasedStaticRegistryExtension;

template<class Extendable, class ...ETags, class ...EDataTypes>
struct TagBasedStaticRegistryExtension<
  Extendable, std::tuple<ETags...>, std::tuple<EDataTypes...>
  >
{
  Extendable & reg_;
  using extension_registry_type = TagBasedStaticRegistry<
    std::tuple<ETags...>, std::tuple<EDataTypes...> >;
  extension_registry_type newReg_;

  template<class ...CArgs>
  TagBasedStaticRegistryExtension(Extendable & reg, CArgs && ... cargs)
    : reg_(reg), newReg_(std::forward<CArgs>(cargs)...){}

  template<class TagToFind>
  static constexpr bool contains(){
    return Extendable::template contains<TagToFind>() ||
      extension_registry_type::template contains<TagToFind>();
  }

  template<class TagToFind, class T>
  void set(T && o){
    if constexpr(Extendable::template contains<TagToFind>()){
      reg_.template set<TagToFind>(std::forward<T>(o));
    }
    else{
      newReg_.template set<TagToFind>(std::forward<T>(o));
    }
  }

  template<class TagToFind>
  auto & get(){
    if constexpr(Extendable::template contains<TagToFind>()){
      return reg_.template get<TagToFind>();
    }
    else{
      return newReg_.template get<TagToFind>();
    }
  }

  template<class TagToFind>
  const auto & get() const {
    if constexpr(Extendable::template contains<TagToFind>()){
      return reg_.template get<TagToFind>();
    }
    else{
      return newReg_.template get<TagToFind>();
    }
  }
};

template<
  class Tag, class DataType,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<Tag>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<Tag>, std::tuple<DataType>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2,
  class D1, class D2,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2>, std::tuple<D1,D2>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2, class T3,
  class D1, class D2, class D3,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T3>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2,T3>, std::tuple<D1,D2,D3>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2, class T3, class T4,
  class D1, class D2, class D3, class D4,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T3>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T4>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2,T3,T4>, std::tuple<D1,D2,D3,D4>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

// this tag is inside the impl namespace because we do NOT want
// to expose it outside, this is an impl detail
struct SystemTag{};

// this represent ths Q^T *r for QR gauss newton
struct QTransposeResidualTag{};


enum class InternalDiagnostic{
  correctionAbsoluteRelativel2Norm,
  residualAbsoluteRelativel2Norm,
  gradientAbsoluteRelativel2Norm,
  objectiveAbsoluteRelative,
  invalid
};

bool is_absolute_diagnostic(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm: return false;
    case Diagnostic::correctionAbsolutel2Norm: return true;
    case Diagnostic::residualRelativel2Norm:   return false;
    case Diagnostic::residualAbsolutel2Norm:   return true;
    case Diagnostic::gradientRelativel2Norm:   return false;
    case Diagnostic::gradientAbsolutel2Norm:   return true;
    case Diagnostic::objectiveAbsolute:   return true;
    case Diagnostic::objectiveRelative:   return false;
    default: return true;
    };
};

std::string diagnostic_to_log_symbol(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm: return "||delta||(r)";
    case Diagnostic::correctionAbsolutel2Norm: return "||delta||(a)";
    case Diagnostic::residualRelativel2Norm:   return "||R||(r)";
    case Diagnostic::residualAbsolutel2Norm:   return "||R||(a)";
    case Diagnostic::gradientRelativel2Norm:   return "||g||(r)";
    case Diagnostic::gradientAbsolutel2Norm:   return "||g||(a)";
    case Diagnostic::objectiveAbsolute:   return "obj(a)";
    case Diagnostic::objectiveRelative:   return "obj(r)";
    default:  return "invalid";
    };
};

std::string diagnostic_to_log_format(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm:
    case Diagnostic::correctionAbsolutel2Norm:
    case Diagnostic::residualRelativel2Norm:
    case Diagnostic::residualAbsolutel2Norm:
    case Diagnostic::gradientRelativel2Norm:
    case Diagnostic::gradientAbsolutel2Norm:
    case Diagnostic::objectiveAbsolute:
    case Diagnostic::objectiveRelative:
      return "{:.6e}";

    default:
      return "invalid";
    };
};

std::string diagnostic_to_string(InternalDiagnostic d){
  switch(d)
    {
    case InternalDiagnostic::correctionAbsoluteRelativel2Norm:
      return "correctionAbsoluteRelativel2Norm";
    case InternalDiagnostic::residualAbsoluteRelativel2Norm:
      return "residualAbsoluteRelativel2Norm";
    case InternalDiagnostic::gradientAbsoluteRelativel2Norm:
      return "gradientAbsoluteRelativel2Norm";
    case InternalDiagnostic::objectiveAbsoluteRelative:
      return "objectiveAbsoluteRelative";
    default:
      return "invalid";
    };
};

InternalDiagnostic public_to_internal_diagnostic(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm:
    case Diagnostic::correctionAbsolutel2Norm:
      return InternalDiagnostic::correctionAbsoluteRelativel2Norm;

    case Diagnostic::residualRelativel2Norm:
    case Diagnostic::residualAbsolutel2Norm:
      return InternalDiagnostic::residualAbsoluteRelativel2Norm;

    case Diagnostic::gradientRelativel2Norm:
    case Diagnostic::gradientAbsolutel2Norm:
      return InternalDiagnostic::gradientAbsoluteRelativel2Norm;

    case Diagnostic::objectiveRelative:
    case Diagnostic::objectiveAbsolute:
      return InternalDiagnostic::objectiveAbsoluteRelative;

    default:
      return InternalDiagnostic::invalid;
    };
}

Diagnostic stop_criterion_to_public_diagnostic(const Stop & sc)
{
  switch (sc)
    {
    case Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance: return Diagnostic::correctionAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfCorrectionBelowTolerance: return Diagnostic::correctionRelativel2Norm;
    case Stop::WhenAbsolutel2NormOfResidualBelowTolerance: return Diagnostic::residualAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfResidualBelowTolerance: return Diagnostic::residualRelativel2Norm;
    case Stop::WhenAbsolutel2NormOfGradientBelowTolerance: return Diagnostic::gradientAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfGradientBelowTolerance: return Diagnostic::gradientRelativel2Norm;
    case Stop::WhenAbsoluteObjectiveBelowTolerance: return Diagnostic::objectiveAbsolute;
    case Stop::WhenRelativeObjectiveBelowTolerance: return Diagnostic::objectiveRelative;
    default: return Diagnostic::invalid;
    };
}

//---------------------------------------------------
//---------------------------------------------------

template<class ValueType>
class InternalDiagnosticDataWithAbsoluteRelativeTracking{
public:
  using value_type = ValueType;

private:
  InternalDiagnostic name_;
  value_type absolute_  = {};
  value_type relative_  = {};
  value_type reference_ = {};

public:
  explicit InternalDiagnosticDataWithAbsoluteRelativeTracking(InternalDiagnostic name)
    : name_(name){}

  InternalDiagnostic name() const{ return name_; }
  value_type getAbsolute() const{ return absolute_; }
  value_type getRelative() const{ return relative_; }

  void update(value_type absoluteValue, bool setAsReferenceValue){
    absolute_ = absoluteValue;
    if (setAsReferenceValue){
      reference_ = absolute_;
      relative_  = absolute_/reference_;
    }
    else{
      relative_ = absolute_/reference_;
    }
  }
};

template<class T>
class DiagnosticsContainer {
  std::vector<T> data_ = {};
  std::vector<Diagnostic> publicNames_ = {};
  std::vector<InternalDiagnostic> internalNames_ = {};
  std::unordered_map<Diagnostic, int> publicNameToIndexMap_ = {};

public:
  explicit DiagnosticsContainer(const std::vector<Diagnostic> & dIn)
    : publicNames_(dIn){ prepare(); }

private:
  void prepare()
  {
    internalNames_.clear();
    data_.clear();
    publicNameToIndexMap_.clear();

    // remove duplicates if needed
    auto last = std::unique(publicNames_.begin(), publicNames_.end());
    publicNames_.erase(last, publicNames_.end());

    // map the public diagnostics to internal ones
    int count = 0;
    for (auto & it : publicNames_){
      const auto internalName = public_to_internal_diagnostic(it);

      // try to see if current value is present
      const auto iter = std::find(internalNames_.begin(), internalNames_.end(), internalName);
      if (iter != internalNames_.end()){
	publicNameToIndexMap_.insert({it, std::distance(internalNames_.begin(), iter)});
      }
      else{
	internalNames_.push_back(internalName);
	data_.push_back( T{internalName} );
	publicNameToIndexMap_.insert({it, count++});
      }
    }
  }

public:
  auto begin()        { return data_.begin(); }
  auto end()          { return data_.end(); }
  auto cbegin() const { return data_.cbegin(); }
  auto cend() const   { return data_.cend(); }

  auto & data() { return data_; }
  const auto & publicNames() const{ return publicNames_; }

  bool contains(const Diagnostic & d) const {
    const auto iter = std::find(publicNames_.begin(), publicNames_.end(), d);
    return iter != publicNames_.end();
  }

  bool contains(const InternalDiagnostic & d) const {
    const auto iter = std::find(internalNames_.begin(), internalNames_.end(), d);
    return iter != internalNames_.end();
  }

  void addIfUnsupported(const Diagnostic & d){
    if (!contains(d)){
      publicNames_.push_back(d);
      prepare();
    }
  }

  int countPublicDiagnostics() const{ return publicNames_.size(); }

  T const & operator[](const Diagnostic & d) const {
    return data_[publicNameToIndexMap_.at(d)];
  }

  T & operator[](const InternalDiagnostic & d){
    const auto iter = std::find(internalNames_.begin(), internalNames_.end(), d);
    const int index = std::distance(internalNames_.begin(), iter);
    return data_[index];
  }

public:
  void print()
  {
    std::cout << "publicNames_\n";
    for (int i=0; i<publicNames_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(publicNames_[i]) << std::endl;
    }
    std::cout << "\n";

    std::cout << "internalNames_\n";
    for (int i=0; i<internalNames_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(internalNames_[i]) << std::endl;
    }
    std::cout << "\n";

    std::cout << "data_\n";
    for (int i=0; i<data_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(data_[i].name()) << std::endl;
    }
    std::cout << "\n";

    std::cout << "dicIndex\n";
    for (auto & it : publicNameToIndexMap_){
      std::cout << diagnostic_to_string(it.first) << " " << it.second << std::endl;
    }
    std::cout << "\n";
  }
};

class DiagnosticsLogger
{
  std::string rootWithLabels_;
  std::string rootWithoutLabels_;

public:
  DiagnosticsLogger() = default;

  DiagnosticsLogger(const std::vector<Diagnostic> & names){
    resetFor(names);
  }

  void resetFor(const std::vector<Diagnostic> & names)
  {
    rootWithoutLabels_ = "{:2d}  ";
    rootWithLabels_ = "nonlinIter = {:2d}: ";
    for (auto & it : names)
      {
	const auto symbol = diagnostic_to_log_symbol(it);
	const auto format = diagnostic_to_log_format(it);
	rootWithLabels_ += symbol + " = " + format + "  ";
	rootWithoutLabels_ += format + "  ";
      }
  }

  template<class T> void operator()(int iStep, const T & dm) const{
    logTrampoline(iStep, dm);
  }

private:
  template<class T>
  void logTrampoline(int iStep, const T & dm) const
  {
    const int count = dm.countPublicDiagnostics();
    if      (count == 1) { logImpl(iStep, dm, std::make_index_sequence<1u>{}); }
    else if (count == 2) { logImpl(iStep, dm, std::make_index_sequence<2u>{}); }
    else if (count == 3) { logImpl(iStep, dm, std::make_index_sequence<3u>{}); }
    else if (count == 4) { logImpl(iStep, dm, std::make_index_sequence<4u>{}); }
    else if (count == 5) { logImpl(iStep, dm, std::make_index_sequence<5u>{}); }
    else if (count == 6) { logImpl(iStep, dm, std::make_index_sequence<6u>{}); }
    else if (count == 7) { logImpl(iStep, dm, std::make_index_sequence<7u>{}); }
    else if (count == 8) { logImpl(iStep, dm, std::make_index_sequence<8u>{}); }
    else if (count == 9) { logImpl(iStep, dm, std::make_index_sequence<9u>{}); }
    else if (count == 10){ logImpl(iStep, dm, std::make_index_sequence<10u>{}); }
    else if (count == 11){ logImpl(iStep, dm, std::make_index_sequence<11u>{}); }
    else{ std::runtime_error("count not implemented yet");}
    ;

  }

  template <class T, std::size_t ... Is>
  void logImpl(int iStep,
	       const DiagnosticsContainer<
	         InternalDiagnosticDataWithAbsoluteRelativeTracking<T>
	       > & dm,
	       std::index_sequence<Is...> const &) const
  {
    using type = InternalDiagnosticDataWithAbsoluteRelativeTracking<T>;
    auto lam = [](Diagnostic d, const type & data) -> typename type::value_type {
        return is_absolute_diagnostic(d) ? data.getAbsolute() : data.getRelative();
    };

    //(..., (std::cout << lam(pd[Is], dm[pd[Is]]) << "\n"));
    const auto pd = dm.publicNames();
    PRESSIOLOG_INFO(rootWithLabels_, iStep, lam(pd[Is], dm[pd[Is]]) ...);
  }
};

//---------------------------------------------------

struct NewtonTag{};
struct GaussNewtonNormalEqTag{};
struct WeightedGaussNewtonNormalEqTag{};
struct LevenbergMarquardtNormalEqTag{};
struct GaussNewtonQrTag{};

// =====================================================
// =====================================================

template<class T>
class LevenbergMarquardtDamping
{
  static_assert(std::is_floating_point<T>::value, "");
  using value_type = T;
  value_type v_{1};
public:
  LevenbergMarquardtDamping & operator = (T v){ v_ = v; return *this; }
  LevenbergMarquardtDamping & operator *= (T v){ v_ *= v; return *this; }
  operator value_type () const { return v_; }
};

template<class RegistryType, class SystemType>
void compute_residual_and_jacobian(RegistryType & reg,
				   const SystemType & system)
{
  const auto & state = reg.template get<StateTag>();
  auto & r = reg.template get<ResidualTag>();
  auto & j = reg.template get<JacobianTag>();
#if defined PRESSIO_ENABLE_CXX20
  if constexpr(SystemWithResidualAndJacobian<SystemType>)
#else
  if constexpr(SystemWithResidualAndJacobian<SystemType>::value)
#endif
  {
    system.residual(state, r);
    system.jacobian(state, j);
  }
  else{
    system.residualAndJacobian(state, r, j);
  }
}

template<class RegistryType>
void compute_gradient(RegistryType & reg)
{
  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();

  const auto & r = reg.template get<ResidualTag>();
  const auto & j = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  // compute gradient (g_ = J_r^T r)
  ::pressio::ops::product(pT, 1, j, r, 0, g);
}

template<class T>
auto compute_half_sum_of_squares(const T & operand)
{
  static_assert(Traits<T>::rank == 1, "");
  const auto normVal = ::pressio::ops::norm2(operand);
  using sc_type = mpl::remove_cvref_t<decltype(normVal)>;
  constexpr auto one  = ::pressio::utils::Constants<sc_type>::one();
  constexpr auto two  = ::pressio::utils::Constants<sc_type>::two();
  return std::pow(normVal, two)*(one/two);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(GaussNewtonNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  auto & r = reg.template get<ResidualTag>();
  system.residual(state, r);
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(GaussNewtonQrTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  auto & r = reg.template get<ResidualTag>();
  system.residual(state, r);
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(LevenbergMarquardtNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  auto & r = reg.template get<ResidualTag>();
  system.residual(state, r);
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(WeightedGaussNewtonNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  const auto & W = reg.template get<WeightingOperatorTag>();
  auto & r  = reg.template get<ResidualTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  system.residual(state, r);
  W.get()(r, Wr);

  const auto v = ::pressio::ops::dot(r, Wr);
  using sc_t = mpl::remove_cvref_t< decltype(v) >;
  constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
  constexpr auto two  = ::pressio::utils::Constants<sc_t>::two();
  return v*(one/two);
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires (   RealValuedSystemWithResidualAndJacobian<SystemType>
	  || RealValuedSystemWithFusedResidualAndJacobian<SystemType>)
#else
template<
  class RegistryType, class SystemType,
  mpl::enable_if_t<
       RealValuedSystemWithResidualAndJacobian<SystemType>::value
    || RealValuedSystemWithFusedResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(GaussNewtonNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{

  compute_residual_and_jacobian(reg, system);

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  auto & H = reg.template get<HessianTag>();
  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  // H = J_r^T J_r, g = J_r^T r
  ::pressio::ops::product(pT, pnT, 1, J, 0, H);
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  return compute_half_sum_of_squares(r);
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires (   RealValuedSystemWithResidualAndJacobian<SystemType>
	  || RealValuedSystemWithFusedResidualAndJacobian<SystemType>)
#else
template<
  class RegistryType, class SystemType,
  mpl::enable_if_t<
       RealValuedSystemWithResidualAndJacobian<SystemType>::value
    || RealValuedSystemWithFusedResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(GaussNewtonQrTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{

  compute_residual_and_jacobian(reg, system);

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();

  auto & g = reg.template get<GradientTag>();
  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  // g = J_r^T r
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  return compute_half_sum_of_squares(r);
}


#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires (   RealValuedSystemWithResidualAndJacobian<SystemType>
	  || RealValuedSystemWithFusedResidualAndJacobian<SystemType>)
#else
template<
  class RegistryType, class SystemType,
  mpl::enable_if_t<
       RealValuedSystemWithResidualAndJacobian<SystemType>::value
    || RealValuedSystemWithFusedResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(WeightedGaussNewtonNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{
  compute_residual_and_jacobian(reg, system);

  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  const auto & W = reg.template get<WeightingOperatorTag>();
  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  auto & WJ = reg.template get<WeightedJacobianTag>();
  auto & g  = reg.template get<GradientTag>();
  auto & H  = reg.template get<HessianTag>();

  W.get()(r, Wr);
  W.get()(J, WJ);
  ::pressio::ops::product(pT, pnT, 1, J, WJ, 0, H);
  ::pressio::ops::product(pT, 1, J, Wr, 0, g);

  using sc_t = scalar_trait_t<typename SystemType::state_type>;
  const auto v = ::pressio::ops::dot(r, Wr);
  constexpr auto one  = ::pressio::utils::Constants<sc_t>::one();
  constexpr auto two  = ::pressio::utils::Constants<sc_t>::two();
  return v*(one/two);
}


#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires (   RealValuedSystemWithResidualAndJacobian<SystemType>
	  || RealValuedSystemWithFusedResidualAndJacobian<SystemType>)
#else
template<
  class RegistryType, class SystemType,
  mpl::enable_if_t<
       RealValuedSystemWithResidualAndJacobian<SystemType>::value
    || RealValuedSystemWithFusedResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(LevenbergMarquardtNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{
  compute_residual_and_jacobian(reg, system);

  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  auto & H = reg.template get<LevenbergMarquardtUndampedHessianTag>();
  auto & scaledH = reg.template get<HessianTag>();
  const auto & damp = reg.template get<LevenbergMarquardtDampingTag>();

  ::pressio::ops::product(pT, pnT, 1, J, 0, H);
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  // compute scaledH = H + mu*diag(H)
  ::pressio::ops::deep_copy(scaledH, H);
  const auto diagH = ::pressio::diag(H);
  auto diaglmH = ::pressio::diag(scaledH);
  ::pressio::ops::update(diaglmH, 1, diagH, damp);

  return compute_half_sum_of_squares(r);
}

template<class RegistryType>
void solve_newton_step(RegistryType & reg)
{
  /* Newton correction solves: J_r delta = - r
     where: delta = x_k+1 - x_k
     so we solve: J_r (-delta) = r and then scale by -1
  */

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & solver = reg.template get<InnerSolverTag>();
  // solve J_r correction = r
  solver.get().solve(J, r, c);
  // scale by -1 for sign convention
  using c_t = mpl::remove_cvref_t<decltype(c)>;
  using scalar_type = typename ::pressio::Traits<c_t>::scalar_type;
  pressio::ops::scale(c, utils::Constants<scalar_type>::negOne() );
}

template<class RegistryType>
void solve_hessian_gradient_linear_system(RegistryType & reg)
{
  const auto & g = reg.template get<GradientTag>();
  const auto & H = reg.template get<HessianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & solver = reg.template get<InnerSolverTag>();
  solver.get().solve(H, g, c);
}

template<class RegistryType>
void compute_correction(GaussNewtonNormalEqTag /*tag*/,
			RegistryType & reg)
{
  /* Gauss-newton with normal eq, we are solving:
       J_r^T*J_r  (x_k+1 - x_k) = - J_r^T*r
     which we can write as:  H c = g
       H = J_r^T*J_r
       g = J_r^T r
       c = x_k - x_k+1

     IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
     we need to rescale c by -1 after the solve
  */
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(WeightedGaussNewtonNormalEqTag /*tag*/,
			RegistryType & reg)
{
  // this is same as regular GN since we solve H delta = g
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(LevenbergMarquardtNormalEqTag /*tag*/,
			RegistryType & reg)
{
  /*
    For LM we are solving: H c = g
    where:
       H = J_r^T*J_r + lambda*diag(J_r^T J_r)
       g = J_r^T r
       c = x_k - x_k+1

    IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
    we need to rescale c by -1 after
  */
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(GaussNewtonQrTag /*tag*/,
			RegistryType & reg)
{
  /*
    see: https://en.wikipedia.org/wiki/Non-linear_least_squares#QR_decomposition
    but careful that here we use J to refer to jacobian wrt r,
    which is the negative of J_f

     IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
     we need to rescale c by -1 after the solve
  */

  const auto & r  = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & QTr = reg.template get<QTransposeResidualTag>();
  auto & solver = reg.template get<InnerSolverTag>();

  auto QTResid = ::pressio::ops::clone(c);
  // factorize J = QR
  solver.get().computeThin(J);

  // compute Q^T r
  solver.get().applyQTranspose(r, QTr);

  // solve Rfactor c = Q^T r
  solver.get().solve(QTr, c);

  // rescale as said above
  ::pressio::ops::scale(c, -1);
}

// =====================================================
// =====================================================
// =====================================================
// =====================================================

template<class Tag, class RegistryType>
void reset_for_new_solve_loop(Tag /*tag*/, RegistryType & reg){
  // no op
}

template<class RegistryType>
void reset_for_new_solve_loop(LevenbergMarquardtNormalEqTag /*tagJ*/,
			      RegistryType & reg)
{
  auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
  damp = 1;
}

// =====================================================
// =====================================================
// =====================================================
// =====================================================

template<class T, class RegistryType>
void compute_norm_internal_diagnostics(const RegistryType & reg,
				       bool isInitial,
				       InternalDiagnosticDataWithAbsoluteRelativeTracking<T> & metric)
{

  switch(metric.name())
  {
    case InternalDiagnostic::residualAbsoluteRelativel2Norm:
      if constexpr(RegistryType::template contains<ResidualTag>()){
	const auto & operand = reg.template get<ResidualTag>();
	const auto value = ops::norm2(operand);
	metric.update(value, isInitial);
      }
    break;

    case InternalDiagnostic::correctionAbsoluteRelativel2Norm:
      if constexpr(RegistryType::template contains<CorrectionTag>()){
	const auto & operand = reg.template get<CorrectionTag>();
	const auto value = ops::norm2(operand);
	metric.update(value, isInitial);
      }
    break;

    case InternalDiagnostic::gradientAbsoluteRelativel2Norm:
      if constexpr(RegistryType::template contains<GradientTag>()){
	const auto & operand = reg.template get<GradientTag>();
	const auto value = ops::norm2(operand);
	metric.update(value, isInitial);
      }
    break;

    default: return;
  };//end switch
}

//---------------------------------------------------
//---------------------------------------------------
//---------------------------------------------------

class DefaultUpdater{
public:
  template<class RegType, class ... Args>
  void operator()(RegType & reg, Args && ... /*unused*/)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: default update");
    auto & state = reg.template get<StateTag>();
    const auto & c = reg.template get<CorrectionTag>();
    ::pressio::ops::update(state, 1, c, 1);
  }
};

class ArmijoUpdater
{
public:
  template<class RegistryType, class ObjF, class ScalarType>
  void operator()(RegistryType & reg,
		  ObjF objective,
		  ScalarType objectiveValueAtCurrentNewtonStep)
  {
    using scalar_type = std::remove_const_t<ScalarType>;
    PRESSIOLOG_DEBUG("Armijo update");

    // https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf

    // Armijo rule says to backtrack alpha until:
    // f(x_trial) - f(x_k) <= alpha_l * beta * dot(g_k, p_k)
    //
    // where:
    // k = the GN step
    // l = indexes the Armijo backtracking stages
    // p_k is the correction at GN k-th step
    // g_k is the gradient at GN k-th step
    // x_trial = x_k + alpha_l * p_k

    const auto & p_k   = reg.template get<CorrectionTag>();
    const auto & g_k   = reg.template get<GradientTag>();
    auto & state = reg.template get<StateTag>();
    auto & x_trial  = reg.template get<LineSearchTrialStateTag>();

    auto alpha = static_cast<scalar_type>(1);
    constexpr auto alpha_lower_bound = static_cast<scalar_type>(0.001);
    const auto beta = static_cast<scalar_type>(0.0001);
    const auto gkDotpk = ::pressio::ops::dot(g_k, p_k);

    PRESSIOLOG_DEBUG("start backtracking");
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    scalar_type ftrial = {};
    while (true)
      {
	if (std::abs(alpha) <= alpha_lower_bound){
	  const auto msg = ": in armijo, alpha = " + std::to_string(alpha)
	    + " <= " + std::to_string(alpha_lower_bound)
	    + ", too small, exiting line search";
	  throw ::pressio::eh::LineSearchStepTooSmall(msg);
	}

	// update : x_trial = x_k + alpha*p_k
	::pressio::ops::update(x_trial, zero, state, one, p_k, alpha);

	// compute rhs_l = alpha_l * beta * dot(g_k, p_k)
	const auto rhs = alpha * beta * gkDotpk;
	const auto lhs = objective(x_trial) - objectiveValueAtCurrentNewtonStep;
	PRESSIOLOG_DEBUG("alpha = {:5f}: (f_trial-f) = {:6e}, rhs = {:6e}", alpha, lhs, rhs);
	if (lhs <= rhs){
	  PRESSIOLOG_DEBUG("condition satisfied: f_trial-f <= rhs, exiting");
	  // solution update: state = state + alpha*p_k
	  ::pressio::ops::update(state, one, p_k, alpha);
	  break;
	}

	// exit when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
	// change later with some machine epsilon
	if (std::abs(lhs) <= 1e-14){
	  const auto msg = ": in armijo, exiting line search";
	  throw ::pressio::eh::LineSearchObjFunctionChangeTooSmall(msg);
	}

	/* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
	alpha *= 0.5;
      }
  }
};

class BacktrackStrictlyDecreasingObjectiveUpdater
{
public:
  template<class RegistryType, class ObjF, class ScalarType>
  void operator()(RegistryType & reg,
		  ObjF objective,
		  ScalarType objectiveValueAtCurrentNewtonStep)
  {
    using scalar_type = std::remove_const_t<ScalarType>;
    PRESSIOLOG_DEBUG("BacktrackStrictlyDecreasingObjective update");

    auto & trialState  = reg.template get<LineSearchTrialStateTag>();
    auto & state = reg.template get<StateTag>();
    const auto & p_k   = reg.template get<CorrectionTag>();

    auto alpha = static_cast<scalar_type>(1);
    constexpr auto alpha_lower_bound = static_cast<scalar_type>(0.001);
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();

    PRESSIOLOG_DEBUG("start backtracking");
    while (true)
      {
	if (std::abs(alpha) <= alpha_lower_bound){
	  /*
	    Presently set an exit alpha drops below 0.001; anything smaller
	    than this is probably unreasonable. Note that this quantity doesn't depend on
	    the dimension or magnitude of the state
	  */
	  const auto msg = ": in BacktrackStrictlyDecreasingObjective: alpha = "
	    + std::to_string(alpha)
	    + " <= " + std::to_string(alpha_lower_bound)
	    + ", too small, exiting line search";
	  throw ::pressio::eh::LineSearchStepTooSmall(msg);
	}

	// update : trialState = state + alpha*p_k
	::pressio::ops::update(trialState, zero, state, one, p_k, alpha);
	if (objective(trialState) <= objectiveValueAtCurrentNewtonStep){
	  PRESSIOLOG_DEBUG("backtrack; condition satisfied with alpha= {:6e}", alpha);
	  ::pressio::ops::update(state, one, p_k, alpha);
	  break;
	}

	/* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
	alpha *= 0.5;
      }
  }
};

template<class RegistryType, class ObjF, class ScalarType, class StateType>
auto lm_gain_factor(RegistryType & reg,
		    ObjF & objective,
		    ScalarType objectiveValueAtCurrentNewtonStep,
		    StateType & cDiagH)
{
  using scalar_type = std::remove_const_t<ScalarType>;
  constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
  constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
  constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();

  const auto & state = reg.template get<StateTag>();
  const auto & correction  = reg.template get<CorrectionTag>();
  auto & trialState  = reg.template get<LineSearchTrialStateTag>();
  const auto & g = reg.template get<GradientTag>();
  const auto & H = reg.template get<LevenbergMarquardtUndampedHessianTag>();
  const auto & damp = reg.template get<LevenbergMarquardtDampingTag>();

  // numerator
  ::pressio::ops::update(trialState, zero, state, one, correction, one);
  const auto num = objective(trialState) - objectiveValueAtCurrentNewtonStep;

  // denominator
  const auto diagH = ::pressio::diag(H);
  ::pressio::ops::elementwise_multiply(one, correction, diagH, zero, cDiagH);
  const auto den1 = ::pressio::ops::dot(correction, cDiagH);
  const auto den2 = ::pressio::ops::dot(correction, g);
  const auto denom = (one/two)*(damp*den1 + den2);

  return num/denom;
}


template<class ScalarType, class StateType>
class LMSchedule1Updater
{
  using scalar_type = std::remove_const_t<ScalarType>;
  using cnst = pressio::utils::Constants<scalar_type>;
  const scalar_type beta_     = cnst::two();
  const scalar_type gammaInv_ = cnst::one()/cnst::three();
  const scalar_type p_ = cnst::three();
  const scalar_type tau_ = cnst::one();
  scalar_type nu_ = cnst::two();
  StateType cDiagH_;

public:
  LMSchedule1Updater(const StateType & stateIn) : cDiagH_(ops::clone(stateIn)){}

public:
  template<class RegType, class Objective>
  void operator()(RegType & reg,
		  Objective obj,
		  scalar_type objectiveValueAtCurrentNewtonStep)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm1 update");

    auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
    const auto & correction  = reg.template get<CorrectionTag>();
    auto & state = reg.template get<StateTag>();

    const scalar_type rho = lm_gain_factor(reg, obj, objectiveValueAtCurrentNewtonStep, cDiagH_);
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();
    if (rho > 0){
      PRESSIOLOG_DEBUG("lm1 update: rho>0");
      ::pressio::ops::update(state, one, correction, one);
      const scalar_type a = one - (beta_ - one)*std::pow(two*rho - one, p_);
      damp *= std::max(gammaInv_, a);
      nu_ = beta_;
      PRESSIOLOG_DEBUG("lm1 update: rho>0, {}, {}", nu_, damp);
    }
    else{
      damp *= nu_;
      nu_ *= two;
      PRESSIOLOG_DEBUG("lm1 update: rho<=0, {}, {}", nu_, damp);
    }
  }
};

template<class ScalarType, class StateType>
class LMSchedule2Updater
{
  using scalar_type = std::remove_const_t<ScalarType>;

  using cnst		   = pressio::utils::Constants<scalar_type>;
  const scalar_type rho1_	   = static_cast<scalar_type>(0.2);
  const scalar_type rho2_     = static_cast<scalar_type>(0.8);
  const scalar_type beta_	   = cnst::two();
  const scalar_type gammaInv_ = cnst::one()/cnst::three();
  const scalar_type tau_	   = cnst::one();
  StateType cDiagH_;

public:
  LMSchedule2Updater(const StateType & stateIn) : cDiagH_(ops::clone(stateIn)){}

public:
  template<class RegType, class Objective>
  void operator()(RegType & reg,
		  Objective obj,
		  scalar_type objectiveValueAtCurrentNewtonStep)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm2 update");
    auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
    const auto & correction  = reg.template get<CorrectionTag>();
    auto & state = reg.template get<StateTag>();

    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto ten  = static_cast<scalar_type>(10);
    constexpr auto seven  = static_cast<scalar_type>(7);
    constexpr auto negSeven  = ::pressio::utils::Constants<scalar_type>::negOne()*seven;
    const auto tenToSev  = std::pow(ten, seven);
    const auto tenToNegSev  = std::pow(ten, negSeven);

    const scalar_type rho = lm_gain_factor(reg, obj, objectiveValueAtCurrentNewtonStep, cDiagH_);
    if (rho < rho1_){
      damp = std::min(damp*beta_, tenToSev);
    }

    if (rho > rho2_){
      damp = std::max( tenToNegSev, damp*gammaInv_);
    }

    if (rho > 0){
      ::pressio::ops::update(state, one, correction, one);
    }
  }
};


//---------------------------------------------------
//---------------------------------------------------

template<
  class ProblemTag,
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class NormDiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void root_solving_loop_impl(ProblemTag /*problemTag*/,
			    const UserDefinedSystemType & system,
			    RegistryType & reg,
			    Stop stopEnumValue,
			    ToleranceType stopTolerance,
			    NormDiagnosticsContainerType & normDiagnostics,
			    const DiagnosticsLoggerType & logger,
			    int maxIters,
			    UpdaterType && updater)
{

  using state_type = typename UserDefinedSystemType::state_type;
  auto objective = [&reg, &system](const state_type & stateIn){
    auto & r = reg.template get<ResidualTag>();
    system.residual(stateIn, r);
    return ::pressio::ops::norm2(r);
  };

  auto mustStop = [&normDiag = std::as_const(normDiagnostics),
		   stopEnumValue, maxIters, stopTolerance](int stepCount){
    const Diagnostic stopDiag = stop_criterion_to_public_diagnostic(stopEnumValue);
    switch (stopEnumValue){
      case Stop::AfterMaxIters:
	return stepCount == maxIters;
      default:
	if (is_absolute_diagnostic(stopDiag)){
	  return normDiag[stopDiag].getAbsolute() < stopTolerance;
	}else{
	  return normDiag[stopDiag].getRelative() < stopTolerance;
	}
    };
  };

  int iStep = 0;
  while (++iStep <= maxIters){
    /* stage 1 */
    try{
      compute_residual_and_jacobian(reg, system);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    /* stage 2 */
    solve_newton_step(reg);

    /* stage 3 */
    std::for_each(normDiagnostics.begin(), normDiagnostics.end(),
		  [&reg, iStep](auto & v){
		    const bool isFirstIteration = iStep==1;
		    compute_norm_internal_diagnostics(reg, isFirstIteration, v);
		  });
    logger(iStep, normDiagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    /* stage 5 */
    try{
      const auto currentObjValue =
	normDiagnostics[InternalDiagnostic::residualAbsoluteRelativel2Norm].getAbsolute();
      updater(reg, objective, currentObjValue);
    }
    catch (::pressio::eh::LineSearchStepTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
    catch (::pressio::eh::LineSearchObjFunctionChangeTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
  }
}


template<class Tag, class StateType, class RegistryType, class NormValueType>
class RootFinder
{
  Tag tag_;
  int maxIters_ = 100;
  Stop stopEnValue_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  NormValueType stopTolerance_ = 0.000001;
  Update updateEnValue_ = Update::Standard;
  RegistryType reg_;
  using norm_diagnostics_container = DiagnosticsContainer<
    InternalDiagnosticDataWithAbsoluteRelativeTracking<NormValueType> >;
  norm_diagnostics_container normDiagnostics_;
  DiagnosticsLogger diagnosticsLogger_ = {};

public:
  RootFinder(Tag tagIn,
	     RegistryType && reg,
	     const std::vector<Diagnostic> & diags)
    : tag_(tagIn), reg_(std::move(reg)), normDiagnostics_(diags)
  {

    // currently we don't have the diagonostics stuff all flushed out
    // so we limit it to work for a specific case
    const auto & publicDiags = normDiagnostics_.publicNames();
    assert(publicDiags.size() == 4);
    assert(publicDiags[0] == Diagnostic::residualAbsolutel2Norm);
    assert(publicDiags[1] == Diagnostic::residualRelativel2Norm);
    assert(publicDiags[2] == Diagnostic::correctionAbsolutel2Norm);
    assert(publicDiags[3] == Diagnostic::correctionRelativel2Norm);

    // check the stop criterion uses a metric already supported in the diagonstics
    const auto stopMetric = stop_criterion_to_public_diagnostic(stopEnValue_);
    normDiagnostics_.addIfUnsupported(stopMetric);
    // need to reset the logger since the names might have changed
    diagnosticsLogger_.resetFor(publicDiags);
  }

  // query/set update criterion
  Update currentUpdateCriterion() const   { return updateEnValue_; }
  void setUpdateCriterion(Update value) { updateEnValue_ = value; }

  // query/set stop criterion, tolerance
  Stop currentStopCriterion() const          { return stopEnValue_; }
  void setStopCriterion(Stop value)	     { stopEnValue_ = value; }
  void setStopTolerance(NormValueType value) { stopTolerance_ = value; }
  void setMaxIterations(int newMax)          { maxIters_ = newMax; }

  template<class SystemType>
  void solve(const SystemType & system, StateType & solutionInOut)
  {
    // deep copy the initial guess
    ::pressio::ops::deep_copy(reg_.template get<InitialGuessTag>(), solutionInOut);

    if (updateEnValue_ == Update::Standard){
      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, StateType &>(reg_, solutionInOut);
      root_solving_loop_impl(tag_, system, extReg,
			     stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_,
			     maxIters_,
			     DefaultUpdater());
    }
    else if (updateEnValue_ == Update::BacktrackStrictlyDecreasingObjective){

      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, LineSearchTrialStateTag,
	StateType &, StateType>(reg_, solutionInOut, system.createState());

      root_solving_loop_impl(tag_, system, extReg,
			     stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_,
			     maxIters_,
			     BacktrackStrictlyDecreasingObjectiveUpdater{});
    }
    else{
      throw std::runtime_error("Invalid criterion");
    }
  }

  // this method can be used when the solver is applied
  // to the same system used for constructing it
  void solve(StateType & solutionInOut)
  {
    auto * system = reg_.template get<SystemTag>();
    assert(system != nullptr);
    this->solve(*system, solutionInOut);
  }
};

template<
  class ProblemTag,
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class DiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void nonlin_ls_solving_loop_impl(ProblemTag problemTag,
				 const UserDefinedSystemType & system,
				 RegistryType & reg,
				 Stop stopEnumValue,
				 ToleranceType stopTolerance,
				 DiagnosticsContainerType & diagnostics,
				 const DiagnosticsLoggerType & logger,
				 int maxIters,
				 UpdaterType && updater)
{

  auto mustStop = [&normDiag = std::as_const(diagnostics),
		   stopEnumValue, maxIters, stopTolerance](int stepCount){
    const Diagnostic stopDiag = stop_criterion_to_public_diagnostic(stopEnumValue);
    switch (stopEnumValue){
    case Stop::AfterMaxIters:
      return stepCount == maxIters;
    default:
      if (is_absolute_diagnostic(stopDiag)){
	return normDiag[stopDiag].getAbsolute() < stopTolerance;
      }else{
	return normDiag[stopDiag].getRelative() < stopTolerance;
      }
    };
  };

  int iStep = 0;
  while (++iStep <= maxIters){
    const bool isFirstIteration = iStep==1;

    // 1. compute operators
    try{
      auto objValue = compute_nonlinearls_operators_and_objective(problemTag, reg, system);
      diagnostics[InternalDiagnostic::objectiveAbsoluteRelative].update(objValue, isFirstIteration);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    // 2. solve for correction
    compute_correction(problemTag, reg);

    /* stage 3 */
    std::for_each(diagnostics.begin(), diagnostics.end(),
		  [&reg, isFirstIteration](auto & v){
		    compute_norm_internal_diagnostics(reg, isFirstIteration, v);
		  });
    logger(iStep, diagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    // 5. run update and continue
    try{
      using state_type = typename UserDefinedSystemType::state_type;
      auto objective = [&reg, &system, problemTag](const state_type & stateIn){
	return compute_nonlinearls_objective(problemTag, reg, stateIn, system);
      };
      const auto currObjValue = diagnostics[InternalDiagnostic::objectiveAbsoluteRelative].getAbsolute();
      updater(reg, objective, currObjValue);

    }
    catch (::pressio::eh::LineSearchStepTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
    catch (::pressio::eh::LineSearchObjFunctionChangeTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
  }
}


template<class Tag, class StateType, class RegistryType, class ScalarType>
class NonLinLeastSquares
{
  static_assert(std::is_floating_point_v<ScalarType>,
		"Impl currently only supporting floating point");

  Tag tag_;
  int maxIters_ = 100;
  Stop stopEnValue_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  ScalarType stopTolerance_ = static_cast<ScalarType>(0.000001);
  Update updateEnValue_ = Update::Standard;

  RegistryType reg_;

  using diagnostics_container = DiagnosticsContainer<
    InternalDiagnosticDataWithAbsoluteRelativeTracking<ScalarType> >;
  diagnostics_container diagnostics_;
  DiagnosticsLogger diagnosticsLogger_ = {};

public:
  NonLinLeastSquares(Tag tagIn,
		     RegistryType && reg,
		     const std::vector<Diagnostic> & diags)
    : tag_(tagIn), reg_(std::move(reg)), diagnostics_(diags)
  {

    // currently we don't have the diagonostics stuff all flushed out
    // so we limit it to work for a specific case
    const auto & publicDiags = diagnostics_.publicNames();
    assert(publicDiags.size() == 6);
    assert(publicDiags[0] == Diagnostic::objectiveAbsolute);
    assert(publicDiags[1] == Diagnostic::objectiveRelative);
    assert(publicDiags[2] == Diagnostic::correctionAbsolutel2Norm);
    assert(publicDiags[3] == Diagnostic::correctionRelativel2Norm);
    assert(publicDiags[4] == Diagnostic::gradientAbsolutel2Norm);
    assert(publicDiags[5] == Diagnostic::gradientRelativel2Norm);

    // check that the stopping criterion uses a metric already
    // supported in the diagonstics, otherwise add it
    const auto stopMetric = stop_criterion_to_public_diagnostic(stopEnValue_);
    diagnostics_.addIfUnsupported(stopMetric);
    // need to reset the logger since the names might have changed
    diagnosticsLogger_.resetFor(publicDiags);
  }

  // query/set update criterion
  Update currentUpdateCriterion() const  { return updateEnValue_; }
  void setUpdateCriterion(Update value)  { updateEnValue_ = value; }

  // query/set stop criterion, tolerance
  Stop currentStopCriterion() const       { return stopEnValue_; }
  void setStopCriterion(Stop value)	  { stopEnValue_ = value; }
  void setStopTolerance(ScalarType value) { stopTolerance_ = value; }
  void setMaxIterations(int newMax)       { maxIters_ = newMax; }

  // this method can be used when passing a system object
  // that is different but syntactically and semantically equivalent
  // to the one used for constructing the solver
  template<class SystemType>
  void solve(const SystemType & system, StateType & solutionInOut)
  {
    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    // if not applicable, this is a noop
    reset_for_new_solve_loop(tag_, reg_);

    switch (updateEnValue_)
    {
      case Update::Standard:
	solve_with_standard_update_impl(system, solutionInOut);
	break;

      case Update::BacktrackStrictlyDecreasingObjective:
      case Update::Armijo:{
	solve_with_line_search_impl(system, solutionInOut);
	break;
      }

      case Update::LMSchedule1:
      case Update::LMSchedule2:{
	if constexpr(std::is_same_v<Tag, LevenbergMarquardtNormalEqTag>){
	  solve_lm_impl(system, solutionInOut);
	  break;
	}
	else{
	  throw std::runtime_error("invalid Update value");
	}
      }

      default:
	throw std::runtime_error("invalid Update value");
    };
  }

  // this method can be used when the solver is applied
  // to the same system used for constructing it
  void solve(StateType & solutionInOut)
  {
    auto * system = reg_.template get<SystemTag>();
    assert(system != nullptr);
    this->solve(*system, solutionInOut);
  }

private:
  template<class SystemType>
  void solve_with_standard_update_impl(const SystemType & system,
				       StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, StateType &>(reg_, solutionInOut);

    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    reset_for_new_solve_loop(tag_, extReg);

    nonlin_ls_solving_loop_impl(tag_, system, extReg,
				stopEnValue_, stopTolerance_,
				diagnostics_, diagnosticsLogger_,
				maxIters_,
				DefaultUpdater());
  }

  template<class SystemType>
  void solve_with_line_search_impl(const SystemType & system,
				   StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, LineSearchTrialStateTag,
      StateType &, StateType>(reg_, solutionInOut, system.createState());

    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    reset_for_new_solve_loop(tag_, extReg);

    if (updateEnValue_ == Update::BacktrackStrictlyDecreasingObjective){
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  BacktrackStrictlyDecreasingObjectiveUpdater{});
    }else{
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  BacktrackStrictlyDecreasingObjectiveUpdater{});
    }
  }

  template<class SystemType>
  void solve_lm_impl(const SystemType & system,
		     StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, LineSearchTrialStateTag,
      StateType &, StateType>(reg_, solutionInOut, system.createState());

    if (updateEnValue_ == Update::LMSchedule1){
      using up_t = LMSchedule1Updater<ScalarType, StateType>;
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  up_t(solutionInOut));
    }else{
      using up_t = LMSchedule2Updater<ScalarType, StateType>;
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  up_t(solutionInOut));
    }
  }
};

}}}
#endif
