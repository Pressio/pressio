
#ifndef PRESSIO_ROM_EPETRA_SKELETON_HPP_
#define PRESSIO_ROM_EPETRA_SKELETON_HPP_

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

namespace pressio{ namespace rom{ namespace test{

class EpetraSkeleton{
protected:
  using nativeVec	= Epetra_Vector;
  using jacobian_type	= Epetra_CrsMatrix;

/* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using velocity_type	= state_type;
  using dense_matrix_type = Epetra_MultiVector;

public:
  EpetraSkeleton() = default;
  ~EpetraSkeleton() = default;

public:
  state_type const & getInitialState() const;

  void velocity(const state_type & u,
		const scalar_type /* t */,
		velocity_type & f) const;

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const;

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const;

  // computes: A = Jac B where B is a multivector
  dense_matrix_type applyJacobian(const state_type & y,
				  const dense_matrix_type & B,
				  scalar_type t) const;

};//end class

}}} //namespace pressio::rom::test
#endif
