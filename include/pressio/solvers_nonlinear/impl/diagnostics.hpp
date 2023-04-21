
#ifndef SOLVERS_NONLINEAR_IMPL_DIAGNOSTICS_HPP_
#define SOLVERS_NONLINEAR_IMPL_DIAGNOSTICS_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

enum class InternalDiagnostic{
  correctionAbsoluteRelativel2Norm,
  residualAbsoluteRelativel2Norm,
  gradientAbsoluteRelativel2Norm,
  objectiveAbsoluteRelative,
  invalid
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

}}}
#endif
