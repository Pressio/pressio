
#ifndef ROM_MASK_DECORATOR_HPP_
#define ROM_MASK_DECORATOR_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace rompp{ namespace rom{ namespace decorator{

/* overload when decorating a residual policy */
template <typename maskable>
class Masked<
  maskable,
  core::meta::enable_if_t<maskable::isResidualPolicy_>
  > : public maskable{

  using typename maskable::fom_rhs_w_t;
  using maskable::fomRhs_;
  using maskable::yFom_;

  mutable std::shared_ptr<fom_rhs_w_t> maskedRhs_ = {};

public:
  Masked() = delete;
  Masked(const maskable & obj) : maskable(obj){}

  template <typename ... Args>
  Masked(Args && ... args)
    : maskable(std::forward<Args>(args)...){}

  ~Masked() = default;

public:
  template <ode::ImplicitEnum odeMethod,  int n,
	     typename ode_state_t,  typename app_t, typename scalar_t>
    fom_rhs_w_t operator()(const ode_state_t & odeY,
			   const std::array<ode_state_t, n> & oldYs,
			   const app_t & app,
			   scalar_t t,
			   scalar_t dt) const
  {
    maskable::template operator()<
      odeMethod, n>(odeY, fomRhs_, oldYs, app, t, dt);

    if (!maskedRhs_){
      maskedRhs_ = std::make_shared<
	fom_rhs_w_t>( app.applyMask(*fomRhs_.data(), t) );
    }
    else
      app.applyMask(*fomRhs_.data(), *maskedRhs_->data(), t);

    return *maskedRhs_;
  }

  template <ode::ImplicitEnum odeMethod,  int n,
	     typename ode_state_t, typename ode_res_t,
	     typename app_t, typename scalar_t>
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t, n> & oldYs,
  		  const app_t & app,
		  scalar_t t,
		  scalar_t dt) const
  {
    maskable::template operator()<
      odeMethod, n>(odeY, fomRhs_, oldYs, app, t, dt);

    app.applyMask(*fomRhs_.data(), *odeR.data(), t);
  }
};//end class



/* overload when decorating a jacobian policy */
template <typename maskable>
class Masked<
  maskable,
  core::meta::enable_if_t<maskable::isResidualPolicy_ == false>
  > : public maskable{

public:
  using typename maskable::apply_jac_return_t;
  mutable std::shared_ptr<apply_jac_return_t> maskedJJ_ = {};

protected:
  using maskable::JJ_;
  using maskable::yFom_;

public:
  Masked() = delete;
  Masked(const maskable & obj) : maskable(obj){}

  template <typename ... Args>
  Masked(Args && ... args)
    : maskable(std::forward<Args>(args)...){}

  ~Masked() = default;

public:
  template <ode::ImplicitEnum odeMethod, typename ode_state_t,
	     typename app_t, typename scalar_t>
  apply_jac_return_t operator()(const ode_state_t & odeY,
				const app_t & app,
				scalar_t t, scalar_t dt) const
  {
    maskable::template operator()<odeMethod>(odeY, JJ_, app, t, dt);
    if (!maskedJJ_){
      maskedJJ_ = std::make_shared<
	apply_jac_return_t
	>( app.applyMask(*JJ_.data(), t) );
    }
    else
      app.applyMask(*JJ_.data(), *maskedJJ_->data(), t);

    return *maskedJJ_;
  }

  template <ode::ImplicitEnum odeMethod, typename ode_state_t,
	     typename ode_jac_t,  typename app_t, typename scalar_t>
  void operator()(const ode_state_t & odeY, ode_jac_t & odeJJ,
  		  const app_t & app, scalar_t t, scalar_t dt) const
  {
    maskable::template operator()<odeMethod>(odeY, JJ_, app, t, dt);
    app.applyMask(*JJ_.data(), *odeJJ.data(), t);
  }

};//end class

}}} //end namespace rompp::rom::decorator

#endif
