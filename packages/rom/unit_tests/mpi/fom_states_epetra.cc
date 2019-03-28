
#include <gtest/gtest.h>
#include "ROM_LSPG"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

using namespace rompp;
using rom_state_t	= core::Vector<Eigen::VectorXd>;
using matrix_w_t	= core::MultiVector<Epetra_MultiVector>;
using decoder_t		= rom::LinearDecoder<matrix_w_t>;
using fom_state_w_t	= core::Vector<Epetra_Vector>;
using fom_state_rec_t	= rom::FomStateReconstructor<fom_state_w_t,decoder_t>;
using fom_states	= rom::FomStatesData<
  fom_state_w_t, 1, fom_state_rec_t>;

struct mytest : fom_states{
  using base_t = fom_states;

  mytest(const fom_state_w_t & y0Fom, const fom_state_rec_t & recObj)
    : base_t(y0Fom, recObj)
  {
    fom_states MyStates(y0Fom, recObj);
    rom_state_t rY(2);
    rY.putScalar(1.2);
    this->reconstructCurrentFomState(rY);
  }

  void check(){
    for (auto i=0; i<this->yFom_.localSize(); i++)
      EXPECT_DOUBLE_EQ(this->yFom_[i], 2.4);
  }
};

TEST(rom_data, fom_states_data){
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map rowMap(6, 0, comm);

  matrix_w_t A(rowMap, 2);
  A.data()->PutScalar(1.0);
  decoder_t decObj(A);

  fom_state_w_t y0Fom(rowMap);

  fom_state_rec_t recStr(y0Fom, decObj);

  mytest Ob(y0Fom, recStr);
  Ob.check();
}
