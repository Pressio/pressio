namespace pressio{ namespace rom{ namespace wls{ namespace impl{
/* 
This header file contains policies for computing the time local residuals; i.e., 
the residual at a single time step,
*/
//============ For the velocity APIs ===========================
struct local_residual_policy_velocityAPI{
public:
  //BDF2
template <typename fom_type, 
          typename fom_state_type, 
          typename residual_type, 
          typename aux_states_container_type, 
          typename scalar_type>
  void operator()(const ::pressio::ode::implicitmethods::BDF2  & tag, 
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  residual_type & residual,
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t,
                  const  scalar_type & dt, 
                  const int & step) const{
    if (step > 0){
    // u^n - 4./3.*u^{n-1} + 1./3.u^{n-2} - 2./3.*dt*f
    appObj.velocity(*yFOM.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<::pressio::ode::implicitmethods::BDF2>(yFOM,residual,auxStatesContainer,dt);
    }
    if (step == 0){
    // u^n - u^{n-1}  - dt*f
    appObj.velocity(*yFOM.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<::pressio::ode::implicitmethods::Euler>(yFOM,residual,auxStatesContainer,dt);
    } 
  }
  // Implicit Euler
template <typename fom_type, 
          typename fom_state_type, 
          typename residual_type, 
          typename aux_states_container_type, 
          typename scalar_type>
  void operator()(const ::pressio::ode::implicitmethods::Euler & tag, 
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  residual_type & residual,
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t,
                  const scalar_type & dt, 
                  const int & step) const{
    appObj.velocity(*yFOM.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<::pressio::ode::implicitmethods::Euler>(yFOM,residual,auxStatesContainer,dt);
  } 
  // Explicit Euler
template <typename fom_type, 
          typename fom_state_type, 
          typename residual_type, 
          typename aux_states_container_type, 
          typename scalar_type>
  void operator()(const ::pressio::ode::explicitmethods::Euler & tag, 
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  residual_type & residual,
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t,
                  const scalar_type & dt, 
                  const int & step) const{
    using nm1 = ode::nMinusOne;
    auto & odeState_nm1 = auxStatesContainer.get(nm1());
    const auto dtnegOne  = ::pressio::utils::constants::negOne<scalar_type>()*dt; //  -1*dt
    constexpr auto negOne  = ::pressio::utils::constants::negOne<scalar_type>(); //  -1*dt
    constexpr auto  One  = ::pressio::utils::constants::one<scalar_type>(); //  -1*dt
    appObj.velocity(*odeState_nm1.data(),t,*residual.data()); // gives f
    ::pressio::containers::ops::do_update(residual,dtnegOne,yFOM,One,odeState_nm1,negOne);
  } 
};

//============================
}}}}
