
#ifndef SOLVERS_NONLINEAR_IMPL_UPDATER_HPP_
#define SOLVERS_NONLINEAR_IMPL_UPDATER_HPP_

template<typename T>
void apply_updater(BaseUpdater* updater,
		  const void* system_as_void,
		  void * reg_as_void,
		  void * diags_as_void)
{
  using system_t = typename T::system_type;
  using reg_t = typename T::reg_type;
  using diag_t = typename T::diag_type;

  auto* p = static_cast<T*>(updater);
  auto* sys = reinterpret_cast<const system_t*>(system_as_void);
  auto* reg = reinterpret_cast<reg_t*>(reg_as_void);
  auto* diag = reinterpret_cast<diag_t*>(diags_as_void);
  p->get()(*sys, *reg, *diag);
}

template<typename T>
void reset_updater(BaseUpdater* updater)
{
  auto* p = static_cast<T*>(updater);
  p->get().reset();
}

struct BaseUpdater
{
  using apply_function_type = void (*)(BaseUpdater*, const void*, void*, void*);
  using reset_function_type = void (*)(BaseUpdater*);
  apply_function_type applyFnc_;
  reset_function_type resetFnc_;

  virtual ~BaseUpdater() = default;

  template <class SystemType, class RegType, class DiagType>
  void operator()(const SystemType & a, RegType & b, DiagType & c){
    (*applyFnc_)(this, &a, &b, &c);
  }

  void reset(){ (*resetFnc_)(this); }
};

template <class SystemType, class RegType, class DiagType, class FunctorType>
class Updater : public BaseUpdater
{
public:
  using functor_type = mpl::remove_cvref_t<FunctorType>;
  using system_type = SystemType;
  using reg_type = RegType;
  using diag_type = DiagType;

private:
  pressio::utils::InstanceOrReferenceWrapper<FunctorType> F_;

public:
  template<class T> Updater(T && Fin) : F_(std::forward<T>(Fin)){}
  ~Updater() = default;
  functor_type & get(){ return F_.get(); }
};


template<typename SystemType, typename RegType, class DiagType>
std::unique_ptr<impl::BaseUpdater>
create_updater(const SystemType & sys,
	       ::pressio::nonlinearsolvers::Update updateE)
{
  using res_t = std::unique_ptr<BaseUpdater>;

  switch (updateE)
    {
    case Update::Standard:{
      using f_t = DefaultUpdater;
      using u_t = Updater<SystemType, RegType, DiagType, f_t>;
      res_t result = pressio::utils::make_unique<u_t>(f_t{});
      result->applyFnc_ = apply_updater<u_t>;
      result->resetFnc_ = reset_updater<u_t>;
      return result;
    }

    // case Update::BacktrackStrictlyDecreasingObjective:{
    //   using f_t = BacktrackStrictlyDecreasingObjectiveUpdater<StateType>;
    //   f_t F(state);
    //   using u_t = Updater<SystemType, StateType, solver_t, f_t>;
    //   res_t result = pressio::utils::make_unique<u_t>(std::move(F));
    //   result->applyFnc_ = apply_updater<u_t>;
    //   result->resetFnc_ = reset_updater<u_t>;
    //   return result;
    // }

    default:
      throw std::runtime_error("Invalid update enum for NewtonRaphson");
      return nullptr;
    }
}

#endif
