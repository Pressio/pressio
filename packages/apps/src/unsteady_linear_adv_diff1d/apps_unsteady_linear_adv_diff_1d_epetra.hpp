
#ifndef ROMPP_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_
#define ROMPP_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_

#include "../../../CORE_ALL"
#include "../apps_ConfigDefs.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include <cmath>
#include "/home/yshimiz/rompp/sources/rompp/packages/apps/src/steady_linear_adv_diff1d/apps_steady_linear_adv_diff_1d_epetra.hpp"
// #include "../steady_linear_adv_diff1d/apps_steady_linear_adv_diff_1d_epetra.hpp"

namespace rompp{ namespace apps{

    class UnsteadyLinAdvDiff1dEpetra: public SteadyLinAdvDiff1dEpetra{

    public:
      UnsteadyLinAdvDiff1dEpetra(Epetra_MpiComm & comm,
				 std::vector<scalar_type> & mu,
				 std::vector<scalar_type> & domain,
				 std::vector<scalar_type> & bc1D)
	: SteadyLinAdvDiff1dEpetra(comm, mu, domain, bc1D){
	
      }
      
      rcp<nativeVec> totaqlResidual() const;
    };
}}
#endif
