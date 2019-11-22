#include <iostream>
#include <string>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_LSPG_UNSTEADY"
#include "APPS_UNSTEADYBURGERS1D"
#include "utils_eigen.hpp"
#include "wls_apis.hpp"
#include "wls_containers.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_basic_meta.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_jacobian_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/solvers/src/meta/solvers_system_has_all_needed_residual_methods.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_sequential_residual_policy_for_residual_api.hpp"
//#include "/Users/ejparis/pressio_repos/pressio/packages/rom/src/experimental/wls/rom_wls_draft.hpp"

// n is the window size
using namespace std;
namespace pressio{ namespace rom{ namespace wls{


template<typename api_tag, typename fom_type, typename wls_state_type,typename ... rest>
struct Specializer{
  using type = void;
};

//struct DefaultApi{};

template<
  typename fom_type, typename wls_state_type,
  typename residual_pol_type, typename jac_pol_type, typename decoder_type
>
struct Specializer<pressio::rom::wls::WlsSystemDefaultApi<fom_type,wls_state_type,residual_pol_type,jac_pol_type,decoder_type>, fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>
{ 
//  static_assert( ::pressio::rom::meta::is_legitimate_residual_policy_for_wls<residual_pol_t>,
//     "You are trying to use Wls with default api but the residual policy passed \
is not admissible for this: maybe you have the wrong api? blas blas");
  using type = pressio::rom::wls::WlsSystemDefaultApi<fom_type, wls_state_type, residual_pol_type, jac_pol_type,decoder_type>;
};
template<
  typename fom_type, typename wls_state_type,
  typename decoder_type, typename fom_state_container_type, typename jtj_jtr_pol_type
>
struct Specializer<pressio::rom::wls::WlsSystemJTJApi<fom_type,wls_state_type,decoder_type,fom_state_container_type, jtj_jtr_pol_type>, fom_type, wls_state_type, decoder_type, fom_state_container_type,jtj_jtr_pol_type>
{ 
  using type = pressio::rom::wls::WlsSystemJTJApi<fom_type, wls_state_type,decoder_type,fom_state_container_type, jtj_jtr_pol_type>;
};

//struct DefaultApi{};
//struct JTJRApi{};
template<typename api_type,typename fom_type,typename wls_state_type,typename ... rest>
using WlsSystem = typename pressio::rom::wls::Specializer<api_type, fom_type, wls_state_type, rest...>::type;
}}} // end namespace rom::wls


int main(int argc, char *argv[]){
  using fom_t		= pressio::apps::Burgers1dEigenResidualApi;
  using scalar_t	= typename fom_t::scalar_type;
  using fom_native_state_t      = typename fom_t::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using rom_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using decoder_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;
  using fom_native_state_t = pressio::apps::Burgers1dEigenResidualApi::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_native_residual_t = pressio::apps::Burgers1dEigenResidualApi::residual_type; 
  using fom_residual_t             = ::pressio::containers::Vector<fom_native_residual_t>;
  using wls_step_jac_t	= pressio::containers::MultiVector<eig_dyn_mat>;
  // Windowed least-squares types. std::vector will be replaced in the future with something custom; e.g., multivector
  using wls_state_t     = std::vector<rom_state_t>; 
  using wls_jacs_t     = std::vector<wls_step_jac_t>;

  std::string checkStr {"PASSED"};

  //---- Information for WLS
  int fomSize = 20;
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  scalar_t dt = 0.01;
  int romSize = 11;
  constexpr std::size_t numStepsInWindow = 5;
  //-------------------------------
  // app object
  fom_t appObj( mu, fomSize);
  auto t0 = static_cast<scalar_t>(0);
  // read from file the jacobian of the decoder
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::rom::test::eigen::readBasis("basis.txt", romSize, fomSize);
  int numBasis = phi.numVectors();
  if( numBasis != romSize ) return 0;
  // create decoder Obj
  decoder_t decoderObj(phi);
  decoderObj.getReferenceToJacobian();
  //using stepper_stencil = ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using fom_state_container_t = typename pressio::ode::StatesContainer<fom_state_t, 2>; 

  // define WLS state
  wls_state_t  wlsState(numStepsInWindow,rom_state_t(romSize) );
  for (int i=0;i<numStepsInWindow; i++){wlsState[i].putScalar(2.0);}
  // define WLS jacobian
  wls_jacs_t  wlsJacs(numStepsInWindow,wls_step_jac_t(fomSize,romSize) );


  fom_state_t yFOM(fomSize);
  fom_state_t yFOM2(fomSize);
  fom_residual_t fom_resid(fomSize);




  using jtr_t = wls_state_t; //this should be the same as WLS state
  using jtj_t = pressio::containers::MultiVector<eig_dyn_mat>;
  jtj_t jtj(romSize,romSize);
  jtr_t jtr(romSize);

  using jtjr_policy_t = pressio::rom::wls::impl::JTJ_JTR_policy_smart<wls_state_t,fom_t,jtj_t,jtr_t,decoder_t>;
  using wls_api_t   = pressio::rom::wls::WlsSystemJTJApi<fom_t,wls_state_t,decoder_t,fom_state_container_t,jtjr_policy_t>;
  using wls_system_t  = pressio::rom::wls::WlsSystem<wls_api_t, fom_t, wls_state_t,decoder_t, fom_state_container_t, jtjr_policy_t>;
  // construct objects and test
  jtjr_policy_t jtjr_policy;
  wls_system_t wls(appObj,jtjr_policy,decoderObj,yFOM,dt,numStepsInWindow);
  wls.JTJ_JTR(wlsState,jtj,jtr); 


//  using resid_t = wls_state_t;
//  using jacobian_t = pressio::containers::MultiVector<eig_dyn_mat>;
//  resid_t residual(romSize);
//  jacobian_t jacobian(romSize,romSize);
//  using residual_policy_t = pressio::rom::wls::impl::residual_policy_naive<wls_state_t,fom_t,resid_t,decoder_t>;
//  using jacobian_policy_t = pressio::rom::wls::impl::jacobian_policy_naive<wls_state_t,fom_t,jacobian_t>;
//  using wls_naive_api_t   = pressio::rom::wls::WlsSystemDefaultApi<fom_t,wls_state_t,residual_policy_t,jacobian_policy_t,decoder_t>;
//  using wls_naive_system_t  = pressio::rom::wls::WlsSystem<wls_naive_api_t, fom_t, wls_state_t,residual_policy_t,jacobian_policy_t,decoder_t>;
//  //construct objects and test
//  residual_policy_t residual_policy;
//  jacobian_policy_t jacobian_policy;
//  wls_naive_system_t wls_naive(appObj,residual_policy,jacobian_policy,decoderObj,yFOM,dt);
  //==================================


  //This works:
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  const fom_state_reconstr_t fom_state_reconstr(yFOM,decoderObj);




  int step = 4;
  scalar_t time = 1.0;
  appObj.timeDiscreteResidual(step,time,dt,*fom_resid.data(),*yFOM.data(),*yFOM2.data()); //nested dependence templated  
  
  decoder_jac_t phi2(fomSize,romSize);
 
  appObj.applyTimeDiscreteJacobian(step,time,dt,*phi.data(),0,*phi2.data(),*yFOM.data(),*yFOM2.data()); //nested dependence templated  


  cout << fom_resid(4) << endl;

  using wls_container_t = pressio::rom::wls::WlsFomStatesContainer<numStepsInWindow,fom_state_t,fom_state_reconstr_t>;
  ///using fom_container_t = pressio::rom::FomStatesStaticContainer<fom_state_t, 2, fom_state_reconstr_t>;
  using fom_container_t = pressio::ode::StatesContainer<fom_state_t, 2>; // just use this fow now

  fom_container_t fom_container(yFOM);
  cout << fom_container[1][0] << endl;
  fom_state_reconstr(wlsState[0],fom_container[1]);
  //fom_container.reconstructCurrentFomState(wlsState[0]); 
  for (int i=1; i< numStepsInWindow; i++){
    fom_state_reconstr(wlsState[i],fom_container[0]);
    fom_state_reconstr(wlsState[i-1],fom_container[1]);
    appObj.timeDiscreteResidual(i,dt*i,dt,*fom_resid.data(),*fom_container[0].data(),*fom_container[1].data());
  }
//  appObj.template timeDiscreteResidual(step,time,dt,fom_resid,y1.data,y2.data); //nested dependence templated  
  using wls_container_t = pressio::rom::wls::WlsFomStatesContainer<numStepsInWindow,fom_state_t,fom_state_reconstr_t>;
  return 0;
}
