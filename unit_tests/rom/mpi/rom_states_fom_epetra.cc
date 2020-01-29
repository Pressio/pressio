
#include <gtest/gtest.h>
#include "ROM_LSPG_UNSTEADY"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

using namespace pressio;
using rom_state_t	= containers::Vector<Eigen::VectorXd>;
using matrix_w_t	= containers::MultiVector<Epetra_MultiVector>;
using decoder_t		= rom::LinearDecoder<matrix_w_t>;
using fom_state_w_t	= containers::Vector<Epetra_Vector>;
using fom_state_rec_t	= rom::FomStateReconstructor<double, fom_state_w_t,decoder_t>;
using fom_states	= rom::FomStatesStaticContainer<fom_state_w_t, 1, fom_state_rec_t>;

struct mytest
{
  fom_states MyStates;

  mytest() = delete;
  mytest(const fom_state_w_t & y0Fom, const fom_state_rec_t & recObj)
    : MyStates{recObj, y0Fom}
  {
    rom_state_t rY(2);
    rY.putScalar(1.2);
    MyStates.reconstructCurrentFomState(rY);
  }

  void check(){
    const auto & yfR = MyStates.getCRefToCurrentFomState();
    const auto sz = yfR.localSize();
    for (auto i=0; i<sz; i++)
      EXPECT_DOUBLE_EQ( yfR[i], 2.4);
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
