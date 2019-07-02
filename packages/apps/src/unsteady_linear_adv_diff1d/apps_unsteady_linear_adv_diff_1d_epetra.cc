
#include "apps_unsteady_linear_adv_diff_1d_epetra.hpp"

#ifdef HAVE_TRILINOS
namespace pressio { namespace apps {
//------------------------------------------------------------------------
// Unsteady setup: Calls the steady setup call before adding initial con
//------------------------------------------------------------------------
void UnsteadyLinAdvDiff1dEpetra::unsteadySetup(){
  setup();
  U0_ = std::make_shared<nativeVec>(*contigMap_);
  U0_->PutScalar(0);
  for (int i = 0; i<nodesPerProc_; i++){
    auto x = (*x_)[i];
    (*U0_)[i] = std::exp(-(x-1)*(x-1)/(2*0.1*0.1));
  }
}

//------------------------------------------------------------------------
// Returns the initial conditions for the states
//------------------------------------------------------------------------
std::shared_ptr<Epetra_Vector>
UnsteadyLinAdvDiff1dEpetra::getInitialState() const{
  return U0_;
}

//------------------------------------------------------------------------
// Returns the velocitys for the unsteady case: signs
//------------------------------------------------------------------------
void UnsteadyLinAdvDiff1dEpetra::velocity(const state_type & u,
					  velocity_type & rhs,
					  const scalar_type /* t*/) const{
  SteadyLinAdvDiff1dEpetra::velocity(u, rhs);
  for (int i = 0; i<nodesPerProc_; i++)
    rhs[i] = -rhs[i];
}

}} //namespace pressio::apps
#endif
