
#ifndef ROM_TEST_LAPLACE1D_EPETRA_HPP_
#define ROM_TEST_LAPLACE1D_EPETRA_HPP_

#include "CORE_ALL"
#include "ODE_ALL"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

class Laplace1dEpetra{
protected:
  using nativeVec = Epetra_Vector;
  template<typename T> using rcp = std::shared_ptr<T>;
  using jacobian_type	= Epetra_CrsMatrix;

/* these types exposed because need to be detected */
public:
  using scalar_type	= double;
  using state_type	= Epetra_Vector;
  using residual_type	= state_type;

public:
  Laplace1dEpetra(Epetra_MpiComm * comm /* other stuff you need */)
    : comm_(comm){}

  ~Laplace1dEpetra() = default;

public:
  void setup();

protected:
  Epetra_MpiComm * comm_;

};

#endif
