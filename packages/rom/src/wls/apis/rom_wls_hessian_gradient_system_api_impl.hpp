
#ifndef ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
#define ROM_WLS_HESSIAN_GRADIENT_SYSTEM_API_IMPL_HPP_
#include "ROM_UTILS"
#include "../time_schemes/rom_wls_select_timescheme_helper.hpp"
#include "../policies/rom_wls_hessian_and_gradient_sequential_policy.hpp"

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename fom_type,
  typename wls_state_type,
  typename decoder_t,
  typename ode_tag,
  typename hessian_t
  >
class SystemHessianGradientApi{

public:
  // aliases for scalar, state, gradient and hessian needed for solver to detect
  // gradient is the same as state for now
  using scalar_type		= typename fom_type::scalar_type;
  using state_type		= wls_state_type;
  using gradient_type		= wls_state_type;
  using hessian_type		= hessian_t;

  using fom_native_state_t	= typename fom_type::state_type;
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t	= pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // policy type (here policy knows how to compute hessian and gradient)
  using hessian_gradient_pol_t	= ::pressio::rom::wls::HessianGradientSequentialPolicy<fom_type,decoder_t>;

  // information on stencil width, time discrete resiudal, time discrete jacobian, etc.
  using time_stencil_t = ::pressio::rom::wls::timeschemes::timescheme_t<ode_tag, fom_state_t, wls_state_type>;
  static constexpr auto timeStencilSize_ = time_stencil_t::state_stencil_size_;

  // in all cases we need types to be wrappers from containers
  static_assert(::pressio::containers::meta::is_vector_wrapper<wls_state_type>::value,
		"For WLS, the state_type must be a pressio container vector");
  static_assert(::pressio::containers::meta::is_matrix_wrapper<hessian_type>::value,
		"WLS: hessian_type must be a pressio container matrix");

  // currently we limit support to eigen types
//  static_assert(::pressio::containers::meta::is_vector_wrapper_eigen<wls_state_type>::value,
//		"WLS: currently supports eigen only");
//  static_assert(::pressio::containers::meta::is_matrix_wrapper_eigen<hessian_type>::value,
//		"WLS: hessian_type must be a Eigen matrix wrapper");

private:
  const fom_type & appObj_;
  const fom_state_reconstr_t  fomStateReconstructorObj_;
  //size of generalized coordinates
  int romSize_			= {};
  // object knowing the time stencil for doing chuncks of step
  time_stencil_t timeSchemeObj_;

  //number of discrete time instances in a window
  int numStepsInWindow_		= {};

  // policy for evaluating the hessian and gradient
  const hessian_gradient_pol_t hessian_gradient_polObj_;

  int activeWindowIndex_	= {};
  // global step number
  int step_s_			= {};
  scalar_type dt_		= {};
  scalar_type windowStartTime_	= {};

  // cache the size of the wls problem: romSize*numStepsInWindow
  int wlsProblemSize_		= romSize_*numStepsInWindow_;
  size_t wlsStencilSize_           = static_cast<size_t>(romSize_*timeStencilSize_);
  // I keep this here since you pass it to the policy. Originally this was owened
  // by the timeSchemeObj but it should be owened here
  wls_state_type wlsStateIC_{static_cast<int>(wlsStencilSize_)};

public:
  SystemHessianGradientApi(const fom_type & appObj,
			   const fom_state_t & yFOM_IC,
			   const fom_state_t & yFOM_Ref,
			   const decoder_t & decoderObj,
			   const int numStepsInWindow)
    : appObj_(appObj),
      fomStateReconstructorObj_(yFOM_Ref, decoderObj),
      romSize_(decoderObj.getReferenceToJacobian().numVectors()),
      timeSchemeObj_(romSize_, yFOM_IC),
      numStepsInWindow_{numStepsInWindow},
      hessian_gradient_polObj_( appObj, yFOM_IC, numStepsInWindow, timeStencilSize_, decoderObj)
  {
    // Set initial condition based on L^2 projection onto trial space
    // note that wlsStateIC_[-romSize:end] contains nm1, wlsStateIC[-2*romSize:-romSize] contains nm2 entry, etc.
    const auto spanStartIndex = romSize_*(timeStencilSize_-2);
  }

  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    // how do we do this for arbitrary types?
    // FRizzi: typically this hessian is small and dense.
    // I think most libraries have this type with a constructor that takes # of rows and cols.
    // As long as we use pressio wrappers, we should be fine I think.
    hessian_type H(wlsProblemSize_, wlsProblemSize_);
    return H;
  }

  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(wlsProblemSize_);
    return g;
  }

  void computeHessianAndGradient(const state_type	      & wls_state,
                                 hessian_type		      & hessian,
                                 gradient_type		      & gradient,
                                 const pressio::solvers::Norm & normType  = ::pressio::solvers::Norm::L2,
                                 scalar_type		      & rnorm = pressio::utils::constants::zero<scalar_type>()) const
  {
    rnorm = pressio::utils::constants::zero<scalar_type>();
    hessian_gradient_polObj_(timeSchemeObj_,
                             wls_state,
                             wlsStateIC_,
                             hessian,
                             gradient,
                             fomStateReconstructorObj_,
                             dt_,
                             numStepsInWindow_,
                             windowStartTime_,
                             step_s_,
                             rnorm);
  }//end computeHessianAndGradient

  // method to advance one window. We may want to put this into some type of window stepper class
  // if we want to have more complex stepping
  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState,
			solverType & solver,
			const int & windowIndex,
			scalar_type dt)
  {
    dt_			= dt; //set time step
    activeWindowIndex_  = windowIndex; //set window number
    windowStartTime_	= windowIndex*dt_*numStepsInWindow_;  //set starting time
    step_s_		= windowIndex*numStepsInWindow_;  //set step number

    solver.solve(*this, wlsState); //solve system

    // Loop to update the the wlsStateIC vector.
    // If we add multistep explicit methods, need to add things here.
    const int start = std::max(0,this->timeStencilSize_ - 1 - this->numStepsInWindow_);
    for (int i = 0; i < start; i++)
    {
      auto wlsTmpState  = ::pressio::containers::span(wlsStateIC_, i*romSize_,     romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC_, (i+1)*romSize_, romSize_);
      for (int k=0; k < romSize_; k++){
	wlsTmpState[k] = wlsTmpState2[k];
      }
    }

    for (int i = start ; i < timeStencilSize_-1; i++)
    {
      auto wlsTmpState  = ::pressio::containers::span(wlsStateIC_, i*romSize_, romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(wlsState, (numStepsInWindow_ - timeStencilSize_+1+i)*romSize_, romSize_);
      for (int k=0; k < romSize_; k++){
	wlsTmpState[k] = wlsTmpState2[k];
      }
    }
    std::cout << " Window " << windowIndex << " completed " << std::endl;
  }// end advanceOneWindow

  // Function to initialize the coefficients for a given FOM IC and reference state
  // computes the ICs via optimal L^2 projection, phi^T phi xhat = phi^T(x - xRef)
  template <typename linear_solver_t>
  void initializeCoeffs(const decoder_t & decoderObj,const fom_state_t & yFOM_IC,const fom_state_t & yRef )
  {
    wls_state_type wlsStateTmp(romSize_);

    ::pressio::rom::utils::setRomCoefficients<linear_solver_t>(decoderObj.getReferenceToJacobian(),wlsStateTmp,yFOM_IC,yRef,romSize_);

    const auto spanStartIndex = romSize_*(timeStencilSize_-2);
    auto wlsInitialStateNm1 = containers::span(wlsStateIC_, spanStartIndex, romSize_);
    for (int i=0; i< this->romSize_; i++){
      wlsInitialStateNm1[i] = wlsStateTmp[i];
    }
  }

};

}}}}//end namespace pressio::rom::src::wls::impl
#endif
