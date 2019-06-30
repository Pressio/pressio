
#ifndef ROMPP_ROM_EPETRA_SKELETON_HPP_
#define ROMPP_ROM_EPETRA_SKELETON_HPP_

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

namespace rompp{ namespace rom{ namespace test{

class EpetraSkeleton{
protected:
  using nativeVec	= Epetra_Vector;
  using jacobian_type	= Epetra_CrsMatrix;

/* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using residual_type	= state_type;

public:
  EpetraSkeleton() = default;
  ~EpetraSkeleton() = default;

public:
  state_type const & getInitialState() const;

  void residual(const state_type & u,
		residual_type & rhs,
		const scalar_type /* t */) const;

  residual_type residual(const state_type & u,
			 const scalar_type t) const;

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const Epetra_MultiVector & B,
		     Epetra_MultiVector & A,
		     scalar_type t) const;

  // computes: A = Jac B where B is a multivector
  Epetra_MultiVector applyJacobian(const state_type & y,
  				   const Epetra_MultiVector & B,
  				   scalar_type t) const;

};//end class

}}} //namespace rompp::rom::test
#endif
