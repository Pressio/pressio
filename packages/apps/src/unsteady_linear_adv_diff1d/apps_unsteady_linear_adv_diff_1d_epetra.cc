
#include "apps_unsteady_linear_adv_diff_1d_epetra.hpp"

namespace rompp { namespace apps {
    void UnsteadyLinAdvDiff1dEpetra::unsteadySetup(){
      U0_ = std::make_shared<nativeVec>(*contigMap_);
      U0_->PutScalar(0);
      for (int i = 0; i<nodesPerProc_; i++)
	{
	  auto GID = MyGlobalNodes_[i];
	  auto x = (*x_)[i];
	  (*U0_)[i] = std::exp(-(x-1)*(x-1)/(2*0.1*0.1));
	}
    }
    
    std::shared_ptr<Epetra_Vector> UnsteadyLinAdvDiff1dEpetra::getInitialState() const{
      return U0_;
    }
}}
