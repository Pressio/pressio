
#ifndef SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_

#include <array>
#include <typeinfo>
#include "solvers_iterative_base.hpp"

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

template<class T, class = void>
struct _has_gradient_instance : std::false_type{};

template<class T>
struct _has_gradient_instance<
  T,
  mpl::enable_if_t<
    mpl::not_void<
      decltype(std::declval<const T&>().gradientCRef())
      >::value
    >
  >
  : std::true_type{};


template<class T, class = void>
struct _has_residual_instance : std::false_type{};

template<class T>
struct _has_residual_instance<
  T,
  mpl::enable_if_t<
    mpl::not_void<
      decltype(std::declval<T const &>().residualCRef() )
      >::value
    >
  >
  : std::true_type{};


template<class T, class = void>
struct _has_correction_instance : std::false_type{};

template<class T>
struct _has_correction_instance<
  T,
  mpl::enable_if_t<
    mpl::not_void<
      decltype(std::declval<T const &>().correctionCRef() )
      >::value
    >
  >
  : std::true_type{};


// ------------------------------------------------------
//
// CORRECTORS
//
// ------------------------------------------------------
template<class StateType, class LinSolverType>
class RJCorrector{
private:
  StateType correction_ = {};
  ::pressio::utils::InstanceOrReferenceWrapper<LinSolverType> solverObj_;

public:
  template <class lsT>
  RJCorrector(const StateType & stateIn, lsT && solverIn)
    : correction_(::pressio::ops::clone(stateIn)),
      solverObj_(std::forward<lsT>(solverIn)){}

public:
  template <class OperatorsT>
  void compute(const OperatorsT & operators)
  {
    PRESSIOLOG_DEBUG("res/jac correction");
    const auto & r = operators.residualCRef();
    const auto & J = operators.jacobianCRef();
    // solve J correction = r
    solverObj_.get().solve(J, r, correction_);
    // scale by -1 for sign convention
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    pressio::ops::scale(correction_, utils::Constants<scalar_type>::negOne() );
  }

  const StateType & correctionCRef() const{ return correction_; }
};

// ------------------------------------------------------
//
// DIAGNOSTICS FNCS and related
//
// ------------------------------------------------------
enum class InternalDiagnostic{
  correctionAbsoluteRelativel2Norm,
  residualAbsoluteRelativel2Norm,
  gradientAbsoluteRelativel2Norm,
  hessianConditionNumber,
  invalid
};

bool is_absolute_diagnostic(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm:
      return false;
    case Diagnostic::correctionAbsolutel2Norm:
      return true;
    case Diagnostic::residualRelativel2Norm:
      return false;
    case Diagnostic::residualAbsolutel2Norm:
      return true;
    case Diagnostic::gradientRelativel2Norm:
      return false;
    case Diagnostic::gradientAbsolutel2Norm:
      return true;
    case Diagnostic::hessianConditionNumber:
      return true;
    default:
      return true;
    };
};

std::string diagnostic_to_string(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm:
      return "correctionRelativel2Norm";
    case Diagnostic::correctionAbsolutel2Norm:
      return "correctionAbsolutel2Norm";
    case Diagnostic::residualRelativel2Norm:
      return "residualRelativel2Norm";
    case Diagnostic::residualAbsolutel2Norm:
      return "residualAbsolutel2Norm";
    case Diagnostic::gradientRelativel2Norm:
      return "gradientRelativel2Norm";
    case Diagnostic::gradientAbsolutel2Norm:
      return "gradientAbsolutel2Norm";
    case Diagnostic::hessianConditionNumber:
      return "hessianConditionNumber";
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
    case InternalDiagnostic::hessianConditionNumber:
      return "hessianConditionNumber";
    default:
      return "invalid";
    };
};

InternalDiagnostic map_to_stored_diag(Diagnostic d){
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

    case Diagnostic::hessianConditionNumber:
      return InternalDiagnostic::hessianConditionNumber;

    default:
      return InternalDiagnostic::invalid;
    };
}

Diagnostic stopping_criterion_to_diagnostic(const Stop & sc)
{
  switch (sc)
    {
    case Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance:
      return Diagnostic::correctionAbsolutel2Norm;

    case Stop::WhenRelativel2NormOfCorrectionBelowTolerance:
      return Diagnostic::correctionRelativel2Norm;

    case Stop::WhenAbsolutel2NormOfResidualBelowTolerance:
      return Diagnostic::residualAbsolutel2Norm;

    case Stop::WhenRelativel2NormOfResidualBelowTolerance:
      return Diagnostic::residualRelativel2Norm;

    case Stop::WhenAbsolutel2NormOfGradientBelowTolerance:
      return Diagnostic::gradientAbsolutel2Norm;

    case Stop::WhenRelativel2NormOfGradientBelowTolerance:
      return Diagnostic::gradientRelativel2Norm;

    default:
      return Diagnostic::invalid;
    };
}

template<class T>
class Metric{
  InternalDiagnostic name_;
  T absolute_  = {};
  T relative_  = {};
  T reference_ = {};

public:
  Metric() = delete;
  explicit Metric(InternalDiagnostic name) : name_(name){}

  InternalDiagnostic name() const{ return name_; }
  T getAbsolute() const{ return absolute_; }
  T getRelative() const{ return relative_; }

  void storeAbsoluteAndUpdateRelative(T value){
    absolute_ = value;
    relative_ = absolute_/reference_;
  }

  void storeAbsoluteAndSetItAsReference(T value){
    absolute_  = value;
    reference_ = absolute_;
    relative_  = absolute_/reference_;
  }
};

// ------------------------------------------------------
//
// DIAGNOSTICS DATA
//
// ------------------------------------------------------
template<class ValueType>
class DiagnosticData
{
  std::vector<Diagnostic> publicDiags_ = {};
  std::vector<InternalDiagnostic> internalDiags_ = {};
  std::vector<Metric<ValueType>> data_ = {};
  std::unordered_map<Diagnostic, int> dicIndex_ = {};
  std::unordered_map<Diagnostic, int> dicSubentry_ = {};

public:
  DiagnosticData(const std::vector<Diagnostic> & dIn)
    : publicDiags_(dIn)
  {
    prepare();
  }

private:
  void prepare()
  {
    internalDiags_.clear();
    data_.clear();
    dicIndex_.clear();
    dicSubentry_.clear();

    // remove duplicates if needed
    auto last = std::unique(publicDiags_.begin(), publicDiags_.end());
    publicDiags_.erase(last, publicDiags_.end());

    // map the public diagnostics to internal ones
    int count = 0;
    for (auto & it : publicDiags_){
      const auto diagIsAbsolute = is_absolute_diagnostic(it);
      const auto storedDiagName = map_to_stored_diag(it);

      // try to see if current value is present
      const auto iter = std::find(internalDiags_.begin(), internalDiags_.end(), storedDiagName);
      if (iter != internalDiags_.end()){
	// if found
	dicIndex_.insert({it, std::distance(internalDiags_.begin(), iter)});
	const auto subEntry = diagIsAbsolute ? 0 : 1;
	dicSubentry_.insert({it, subEntry});
      }
      else{
	internalDiags_.push_back(storedDiagName);
	data_.push_back(Metric<ValueType>{storedDiagName});
	dicIndex_.insert({it, count++});
	const auto subEntry = diagIsAbsolute ? 0 : 1;
	dicSubentry_.insert({it, subEntry});
      }
    }
  }

public:
  void print()
  {
    std::cout << "publicDiags_\n";
    for (int i=0; i<publicDiags_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(publicDiags_[i]) << std::endl;
    }
    std::cout << "\n";

    std::cout << "internalDiags_\n";
    for (int i=0; i<internalDiags_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(internalDiags_[i]) << std::endl;
    }
    std::cout << "\n";
    std::cout << "data_\n";
    for (int i=0; i<data_.size(); ++i){
      std::cout << i << " " << diagnostic_to_string(data_[i].name()) << std::endl;
    }
    std::cout << "\n";

    std::cout << "dicIndex\n";
    for (auto & it : dicIndex_){
      std::cout << diagnostic_to_string(it.first) << " " << it.second << std::endl;
    }
    std::cout << "\n";
    std::cout << "dicSubentry\n";
    for (auto & it : dicSubentry_){
      std::cout << diagnostic_to_string(it.first) << " " << it.second << std::endl;
    }
    std::cout << "\n";
  }

  void addIfUnsupported(const Diagnostic & newD){
    if (!supports(newD)){
      publicDiags_.push_back(newD);
      prepare();
    }
  }

  bool supports(const Diagnostic & d) const {
    const auto iter = std::find(publicDiags_.begin(), publicDiags_.end(), d);
    return iter != publicDiags_.end();
  }

  bool supports(const InternalDiagnostic & d) const {
    const auto iter = std::find(internalDiags_.begin(), internalDiags_.end(), d);
    return iter != internalDiags_.end();
  }

  int publicDiagsCount()     const{ return publicDiags_.size(); }
  const auto & publicDiags() const{ return publicDiags_; }

  auto & data() { return data_; }

  ValueType operator[](const Diagnostic & d) const {
    const auto index = dicIndex_.at(d);
    const auto subentry = dicSubentry_.at(d);
    return (subentry==0) ? data_[index].getAbsolute() : data_[index].getRelative();
  }
};

// ------------------------------------------------------
//
// DIAGNOSTICS EVALUATIONS
//
// ------------------------------------------------------

template<class T, class OperandType>
mpl::enable_if_t< !_has_correction_instance<OperandType>::value >
computeInternalDiagnosticsOnCorrection(bool isInitial,
				       Metric<T> & metric,
				       const OperandType & operand)
{ /*no op */}

template<class T, class OperandType>
mpl::enable_if_t< _has_correction_instance<OperandType>::value >
computeInternalDiagnosticsOnCorrection(bool isInitial,
				     Metric<T> & metric,
				     const OperandType & operand)
{
  switch(metric.name())
    {
    case InternalDiagnostic::correctionAbsoluteRelativel2Norm:
      if (isInitial){
	metric.storeAbsoluteAndSetItAsReference(ops::norm2(operand.correctionCRef()));
      }else{
	metric.storeAbsoluteAndUpdateRelative(ops::norm2(operand.correctionCRef()));
      }
      break;
    };
}

template<class T, class OperandType>
mpl::enable_if_t< !_has_residual_instance<OperandType>::value >
computeInternalDiagnosticsOnResidual(bool isInitial,
				     Metric<T> & metric,
				     const OperandType & operand)
{ /*no op */}

template<class T, class OperandType>
mpl::enable_if_t< _has_residual_instance<OperandType>::value >
computeInternalDiagnosticsOnResidual(bool isInitial,
				     Metric<T> & metric,
				     const OperandType & operand)
{
  switch(metric.name())
    {
    case InternalDiagnostic::residualAbsoluteRelativel2Norm:
      if (isInitial){
	metric.storeAbsoluteAndSetItAsReference(ops::norm2(operand.residualCRef()));
      }else{
	metric.storeAbsoluteAndUpdateRelative(ops::norm2(operand.residualCRef()));
      }
      break;
    };
}

template<class T, class OperandType>
mpl::enable_if_t< !_has_gradient_instance<OperandType>::value >
computeInternalDiagnosticsOnGradient(bool isInitial,
				     Metric<T> & metric,
				     const OperandType & operand)
{ /*no op */}

template<class T, class OperandType>
mpl::enable_if_t< _has_gradient_instance<OperandType>::value >
computeInternalDiagnosticsOnGradient(bool isInitial,
				     Metric<T> & metric,
				     const OperandType & operand)
{
  switch(metric.name())
    {
    case InternalDiagnostic::gradientAbsoluteRelativel2Norm:
      if (isInitial){
	metric.storeAbsoluteAndSetItAsReference(ops::norm2(operand.gradientCRef()));
      }else{
	metric.storeAbsoluteAndUpdateRelative(ops::norm2(operand.gradientCRef()));
      }
      break;
    };
}

// ------------------------------------------------------
//
// DIAGNOSTICS LOGGER
//
// ------------------------------------------------------
class DiagnosticsLogger
{
  std::string rootWithLabels_;
  std::string rootWithoutLabels_;

public:
  DiagnosticsLogger(const std::vector<Diagnostic> & names){
    resetFor(names);
  }

  void resetFor(const std::vector<Diagnostic> & names)
  {
    rootWithoutLabels_ = "{:2d} ";
    rootWithLabels_ = "nonlinIter = {:2d}: ";
    for (auto & it : names)
    {
      switch(it)
	{
	case Diagnostic::correctionAbsolutel2Norm:
	  rootWithLabels_    += "||delta||(a) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;
	case Diagnostic::correctionRelativel2Norm:
	  rootWithLabels_    += "||delta||(r) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;

	case Diagnostic::residualAbsolutel2Norm:
	  rootWithLabels_    += "||R||(a) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;
	case Diagnostic::residualRelativel2Norm:
	  rootWithLabels_    += "||R||(r) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;

	case Diagnostic::gradientAbsolutel2Norm:
	  rootWithLabels_    += "||g||(a) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;
	case Diagnostic::gradientRelativel2Norm:
	  rootWithLabels_    += "||g||(r) = {:.6e} ";
	  rootWithoutLabels_ += "{:.6e} ";
	  break;
	};
    }
  }

  template<class T>
  void log(int iStep, const T & dm)
  {
    logImpl(iStep, dm);
  }

private:
  template<class T>
  void logImpl(int iStep, const T & dm)
  {
    if (dm.publicDiagsCount() == 1){
      logImpl(iStep, dm, std::make_index_sequence<1u>{});
    }
    else if (dm.publicDiagsCount() == 2){
      logImpl(iStep, dm, std::make_index_sequence<2u>{});
    }
    else if (dm.publicDiagsCount() == 3){
      logImpl(iStep, dm, std::make_index_sequence<3u>{});
    }
    else if (dm.publicDiagsCount() == 4){
      logImpl(iStep, dm, std::make_index_sequence<4u>{});
    }
    else if (dm.publicDiagsCount() == 5){
      logImpl(iStep, dm, std::make_index_sequence<5u>{});
    }
    else if (dm.publicDiagsCount() == 6){
      logImpl(iStep, dm, std::make_index_sequence<6u>{});
    }
  }

  template <class T, std::size_t ... Is>
  void logImpl(int iStep, const T  & dm, std::index_sequence<Is...> const &)
  {
    const auto pd = dm.publicDiags();
    PRESSIOLOG_INFO(rootWithLabels_, iStep, dm[pd[Is]] ...);
  }
};

// ------------------------------------------------------
//
// STOPPER
//
// ------------------------------------------------------
template<class T>
class Stopper
{
  T tolerance_ = 0.000001;
  Stop stoppingE_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  std::size_t maxIters_ = {};

public:
  Stopper(std::size_t maxIters) : maxIters_(maxIters){}

  auto criterion() const { return stoppingE_; }
  auto tolerance() const { return tolerance_; }
  void setTolerance(T newTol)   { tolerance_ = newTol; }
  void setCriterion(Stop value) { stoppingE_ = value; }
  void setCriterion(Stop value, T newTol){
    setCriterion(value);
    tolerance_ = newTol;
  }

  bool isDone(std::size_t iStep,
	      const DiagnosticData<T> & dm) const
  {

    switch (stoppingE_)
      {
      case Stop::AfterMaxIters:
	return iStep == maxIters_;

      default:
	const auto d = stopping_criterion_to_diagnostic(stoppingE_);
	return dm[d] < tolerance_;
      };
  }
};

// ------------------------------------------------------
//
// SOLVER
//
// ------------------------------------------------------
template<class StateType, class OperatorsT, class CorrectorT>
class Solver
  : public IterativeBase<Solver<StateType, OperatorsT, CorrectorT>>
{
  using this_type = Solver<StateType, OperatorsT, CorrectorT>;
  using iterative_base_type = IterativeBase<this_type>;
  friend iterative_base_type;
  using typename iterative_base_type::iteration_count_type;
  using norm_value_type = typename OperatorsT::residual_norm_value_type;

  OperatorsT o_;
  CorrectorT c_;
  DiagnosticData<norm_value_type> dm_;
  DiagnosticsLogger logger_;
  Stopper<norm_value_type> stopper_;
  // Updater updater_;

public:
  Solver(OperatorsT && o,
	 CorrectorT && c,
	 const std::vector<Diagnostic> & diags)
    : o_(std::move(o)),
      c_(std::move(c)),
      dm_(diags),
      logger_(dm_.publicDiags()),
      stopper_{iterative_base_type::maxIters_}
  {

    // check that the stopping criterion uses a metric already
    // supported in the diagonstics, otherwise add it
    const auto stopMetric = stopping_criterion_to_diagnostic(stopper_.criterion());
    dm_.addIfUnsupported(stopMetric);
    //dm_.print();
    // need to reset the logger since the names might have changed
    logger_.resetFor(dm_.publicDiags());
  }

public:
  template<typename SystemType>
  void solve(const SystemType & system,
	     StateType & state)
  {
    iteration_count_type iStep = 0;
    bool recomputeSystemJacobian = true;

    while (++iStep <= iterative_base_type::maxIters_)
    {
      recomputeSystemJacobian = true;
      //(iStep_ == 1) ? true : ((iStep_ % jacobianUpdateFreq_) == 0);

      // 1.
      try{
	o_.compute(system, state, recomputeSystemJacobian);
      }
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::NonlinearSolveFailure();
      }
      catch (::pressio::eh::ResidualHasNans const &e){
	PRESSIOLOG_CRITICAL(e.what());
	throw ::pressio::eh::NonlinearSolveFailure();
      }

      // 2.
      c_.compute(o_);

      // 4.
      const bool isFirstIteration = iStep==1;
      auto & metrics = dm_.data();
      std::for_each(metrics.begin(), metrics.end(),
		    [&](auto & m){
		      computeInternalDiagnosticsOnCorrection(isFirstIteration, m, c_);
		      computeInternalDiagnosticsOnResidual(isFirstIteration, m, o_);
		    });

      // 5.
      logger_.log(iStep, dm_);

      // 6.
      if (stopper_.isDone(iStep, dm_)){
	PRESSIOLOG_DEBUG("nonlinsolver: stopping");
	break;
      }

      const auto & correction = c_.correctionCRef();
      using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
      ::pressio::ops::update(state, one, correction, one);
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_SOLVER_HPP_
