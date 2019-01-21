
#ifndef ROM_PRECONDITIONED_MIXIN_HPP_
#define ROM_PRECONDITIONED_MIXIN_HPP_

namespace rompp{ namespace rom{


template <typename T, typename enable = void>
struct is_residual_policy : std::false_type{};

template <typename T>
struct is_residual_policy<
  T,
  typename std::enable_if<
    !std::is_void<
      typename T::app_res_w_t
      >::value
    >::type
  > : std::true_type{};




/* this is for residual */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  core::meta::enable_if_t<
    is_residual_policy<preconditionable>::value
    >
  >{

  using res_w_type  = typename precondionable::app_res_w_type;
  using scalar_type = typename preconditionable::scalar_type;

protected:
  // default constructor
  Preconditioned() = delete;

  Precondtioned(const Preconditionable & obj)
    : myObj_(obj){}

  ~Preconditioned() = default;

public:

  template <::rompp::ode::ImplicitEnum odeMethod,
	     int numAuxStates,
	     typename ode_state_t,
	     typename app_t>
  res_w_type operator()(const ode_state_t & odeY,
			const std::array<ode_state_t,numAuxStates> & oldYs,
			const app_t & app,
			scalar_type t,
			scalar_type dt) const
  {
    res_w_type A = myObj_(odeY, oldYs, app, t, dt);
    app.applyPreconditioner(odeY, A, t);

  }


  template <::rompp::ode::ImplicitEnum odeMethod,
	     int numAuxStates,
	     typename ode_state_t,
	     typename ode_res_t,
	     typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t,numAuxStates> & oldYs,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const
  {
    myObj_(odeY, odeR, oldYs, app, t, dt);
    app.applyPreconditioner(odeY, odeR, t);
  }


private:
  const preconditionable & myObj_;

};//end class


}} //end namespace rompp::rom

#endif
