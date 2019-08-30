
#ifndef PRESSIO_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_
#define PRESSIO_APPS_UNSTEADY_LINEAR_ADV_DIFF_1D_EPETRA_HPP_

#include "../apps_ConfigDefs.hpp"

#ifdef HAVE_TRILINOS
#include "../../../CONTAINERS_ALL"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include <cmath>
#include "../steady_linear_adv_diff1d/apps_steady_linear_adv_diff_1d_epetra.hpp"

namespace pressio{ namespace apps{
class UnsteadyLinAdvDiff1dEpetra: public SteadyLinAdvDiff1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  using nativeMatrix  = Epetra_CrsMatrix;

public:
  UnsteadyLinAdvDiff1dEpetra(const Epetra_MpiComm & comm,
			     const std::vector<scalar_type> & mu,
			     const std::vector<scalar_type> & domain,
			     const std::vector<scalar_type> & bc1D)
    : SteadyLinAdvDiff1dEpetra(comm, mu, domain, bc1D){}
  ~UnsteadyLinAdvDiff1dEpetra() = default;

public:
  void unsteadySetup();

  rcp<nativeVec> getInitialState() const;

  void velocity(const state_type & u,
		const scalar_type /* t*/,
    velocity_type & rhs) const;

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const{
    Epetra_Vector R(*contigMap_);
    velocity(u, t, R);
    return R;
  }

  void applyJacobian(const state_type & y,
		     const Epetra_MultiVector & B,
		     scalar_type /*t*/,
         Epetra_MultiVector &A) const{
    SteadyLinAdvDiff1dEpetra::applyJacobian(y, B, A);
    A.Scale(-1.0);
  }

  Epetra_MultiVector applyJacobian(const state_type &y,
				   const Epetra_MultiVector &B,
				   scalar_type t) const{
    Epetra_MultiVector C( *contigMap_, B.NumVectors());
    applyJacobian(y, B, t, C);
    return C;
  }

protected:
  mutable rcp<nativeVec> U0_;
};

}}
#endif
#endif
