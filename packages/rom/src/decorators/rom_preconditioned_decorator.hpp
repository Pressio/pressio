
#ifndef ROM_PRECONDITIONED_DECORATOR_HPP_
#define ROM_PRECONDITIONED_DECORATOR_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace rompp{ namespace rom{ namespace decorator{

/* overload when decorating a residual policy */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  ::rompp::mpl::enable_if_t<preconditionable::isResidualPolicy_>
  > : public preconditionable{

  using typename preconditionable::fom_rhs_t;
  using preconditionable::yFom_;

public:
  Preconditioned() = delete;
  Preconditioned(const preconditionable & obj) : preconditionable(obj){}

  template <typename ... Args>
  Preconditioned(Args && ... args)
    : preconditionable(std::forward<Args>(args)...){}

  ~Preconditioned() = default;

public:

  //-------------------------------
  // for unsteady LSPG
  //-------------------------------
  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
  >
  fom_rhs_t operator()(const ode_state_t & odeY,
			 const std::array<ode_state_t, n> & oldYs,
			 const app_t & app,
			 scalar_t t,
			 scalar_t dt) const
  {
    auto result = preconditionable::template operator()<
      odeMethod, n>(odeY, oldYs, app, t, dt);

    app.applyPreconditioner(*yFom_.data(), *result.data(), t);
    return result;
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename ode_state_t,
    typename ode_res_t,
    typename app_t,
    typename scalar_t
    >
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t, n> & oldYs,
  		  const app_t & app,
		  scalar_t t,
		  scalar_t dt) const
  {
    preconditionable::template operator()<
      odeMethod, n>(odeY, odeR, oldYs, app, t, dt);

    app.applyPreconditioner(*yFom_.data(), *odeR.data(), t);
  }


  //-------------------------------
  // for steady LSPG
  //-------------------------------
  template <
      typename lspg_state_t,
      typename fom_t>
  fom_rhs_t operator()(const lspg_state_t  & romY,
                         const fom_t   & app) const
  {
    auto result = preconditionable::template operator()<
      lspg_state_t, fom_t>(romY, app);
    app.applyPreconditioner(*yFom_.data(), *result.data());
    return result;
  }

  template <
      typename lspg_state_t,
      typename lspg_residual_t,
      typename fom_t>
  void operator()(const lspg_state_t  & romY,
                  lspg_residual_t & romR,
                  const fom_t   & app) const
  {
    preconditionable::template operator()<
      lspg_state_t, lspg_residual_t, fom_t>(romY, romR, app);
    app.applyPreconditioner(*yFom_.data(), *romR.data());
  }

};//end class



/* overload when decorating a jacobian policy */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  ::rompp::mpl::enable_if_t<preconditionable::isResidualPolicy_ == false>
  > : public preconditionable{

public:
  using typename preconditionable::apply_jac_return_t;

protected:
  using preconditionable::yFom_;

public:
  Preconditioned() = delete;
  Preconditioned(const preconditionable & obj) : preconditionable(obj){}

  template <typename ... Args>
  Preconditioned(Args && ... args)
    : preconditionable(std::forward<Args>(args)...){}

  ~Preconditioned() = default;

public:

  //-------------------------------
  // for unsteady LSPG
  //-------------------------------
  template <
    ode::ImplicitEnum odeMethod,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
  >
  apply_jac_return_t operator()(const ode_state_t & odeY,
				const app_t & app,
				scalar_t t, scalar_t dt) const
  {
    auto JJ = preconditionable::template operator()<odeMethod>(odeY, app, t, dt);
    app.applyPreconditioner(*yFom_.data(), *JJ.data(), t);
    return JJ;
  }

  template <
      ode::ImplicitEnum odeMethod,
      typename ode_state_t,
      typename ode_jac_t,
      typename app_t,
      typename scalar_t>
  void operator()(const ode_state_t & odeY, ode_jac_t & odeJJ,
  		  const app_t & app, scalar_t t, scalar_t dt) const
  {
    preconditionable::template operator()<odeMethod>(odeY, odeJJ, app, t, dt);
    app.applyPreconditioner(*yFom_.data(), *odeJJ.data(), t);
  }


  //-------------------------------
  // for steady LSPG
  //-------------------------------
  template <
      typename lspg_state_t,
      typename app_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
                                const app_t & app) const
  {
    auto JJ = preconditionable::template operator()<
      lspg_state_t, app_t>(romY, app);
    app.applyPreconditioner(*yFom_.data(), *JJ.data());
    return JJ;
  }

  template <
      typename lspg_state_t,
      typename lspg_jac_t,
      typename app_t>
  void operator()(const lspg_state_t & romY,
                  lspg_jac_t & romJJ,
                  const app_t & app) const
  {
    preconditionable::template operator()<
      lspg_state_t, lspg_jac_t, app_t>(romY, romJJ, app);
    app.applyPreconditioner(*yFom_.data(), *romJJ.data());
  }

};//end class

}}} //end namespace rompp::rom::decorator

#endif
