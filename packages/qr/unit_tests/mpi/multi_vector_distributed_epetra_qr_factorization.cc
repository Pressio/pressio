
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "../../src/qr_epetra_multi_vector.hpp"
#include "Epetra_MpiComm.h"

struct epetraMultiVectorR9C4Fixture
  : public ::testing::Test{

public:
  int rank_;
  Epetra_MpiComm * comm_;
  int numProc_;
  const int localSize_ = 3;
  const int numVectors_ = 4;
  int numGlobalEntries_;
  Epetra_Map * dataMap_;
  Epetra_MultiVector * mv_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
    rank_ = comm_->MyPID();
    numProc_ = comm_->NumProc();
    assert(numProc_ == 3);

    numGlobalEntries_ = numProc_ * localSize_;
    dataMap_ = new Epetra_Map(numGlobalEntries_, 0, *comm_);
    mv_ = new Epetra_MultiVector(*dataMap_, numVectors_);
  }

  virtual void TearDown(){
    delete comm_;
    delete dataMap_;
    delete mv_;
  }
};
//-----------------------------------------------------------


TEST_F(epetraMultiVectorR9C4Fixture,
     EpetraMultiVectorQRFactorization){
  using namespace rompp;

  using nat_t = Epetra_MultiVector;
  using mymvec_t = core::MultiVector<nat_t>;
  STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(mymvec_t);

  // fill the matrix
  mymvec_t A( *dataMap_, numVectors_ );
  if(rank_==0){
    A(0,0) = 3.2; A(0,1) = 1.2;  A(0,2) = 1.;
    A(1,0) = 1.2; A(1,2) = -2.2;
    A(2,1) = 4.0; A(2,3) = -2.;
  }
  if(rank_==1){
    A(0,1) = 4.;
    A(1,2) = -1.; A(1,3) = -4.;
    A(2,0) = 0.2; A(2,1) = 5.;  A(2,2) = 1.;
  }
  if(rank_==2){
    A(0,0) = 1.; A(0,1) = 1.1; A(0,2) = 1.25; A(0,3) = -3.;
    A(1,2) = 1.; A(1,1) = 0.1111; A(1,3) = 6.;
  }
  A.data()->Print(std::cout);

  // do QR
  using R_type = rompp::core::Matrix<Eigen::MatrixXd>;
  rompp::qr::QRSolver<mymvec_t, rompp::core::MultiVector, R_type> qrObj;
  qrObj.compute(A);
  //  using Q_t = rompp::core::MultiVector<Epetra_MultiVector>;
  const auto & Q = qrObj.cRefQFactor();
  const auto & R = qrObj.cRefRFactor();
  Q.data()->Print(std::cout << std::setprecision(4));
  if (rank_==0)
    std::cout << "\n" << *R.data();


  Eigen::MatrixXd trueQ(9,9);
  trueQ << -0.897235446547271, -0.039024431200404,  0.173692309541681, -0.034843851055998,
    0.06844323641866,  -0.162666577780091, -0.355064595330458, -0.069618789726194,  0.,
    ///first row
    -0.336463292455227, -0.074296513246923, -0.7624257912921,   -0.286577009249518,
    -0.15470306042946,   0.249902247740785,  0.353018563741718,  0.082744998721086,  0.,
    //end row
    -0.,                 0.530332013749077, -0.155702484341873,  0.126884420222574,
    -0.42395658892426,  -0.401481903098228, -0.154500472643382,  0.560006903183078,  0.,
    //end row
    -0.,                 0.530332013749077, -0.155702484341873, -0.139879089512256,
    0.326913228407908, -0.468149250960295,  0.260689856536577, -0.530040895116606,  0.,
    //end row
    -0.,                 0.,                -0.32818521819755,   0.421198101926496,
    0.751758650006291,  0.022892491611334, -0.138183613491497,  0.360730891969511,  0.,
    //end row
    -0.056077215409204,  0.650532264978526,  0.126820727560628, -0.069095509456248,
    0.063878165111811,  0.719648521641246, -0.169043632892447, -0.04748965432464,   0.,
    //end row
    -0.280386077046022,  0.083927542741894,  0.333731413505014,  0.469211836369866,
    -0.046150317046723,  0.0767206472791,    0.746393155145893,  0.132984059523448,  0.,
    //end row
    -0.,                 0.014729971681881,  0.323860581694954, -0.691846753372529,
    0.336778411839413, -0.08020531631897,   0.229574011030821,  0.493648258802424,  0.,
    //end row
    -0., 0., 0.,-0., 0., 0., 0., 0.,  1.;

  Eigen::MatrixXd trueR(9,4);
  trueR << -3.566510900025401, -1.665493297653372, -0.563576014862505,  0.841158231138066,//endrow
	    0.,  7.542444914314701,  0.8945995630306  , -1.224066825632553,//endrow
	    0.,  0.               ,  3.047059844718702,  2.566115091128629,//endrow
	    0.,  0.               ,  0.               , -7.497277277495907,//endrow
	    0.,  0.               ,  0.               ,  0.               ,//endrow
	    0.,  0.               ,  0.               ,  0.               ,//endrow
	    0.,  0.               ,  0.               ,  0.               ,//endrow
	    0.,  0.               ,  0.               ,  0.               ,//endrow
	    0.,  0.               ,  0.               ,  0.               ;

  // check R factor
  for (auto i=0; i<9; i++)
    for (auto j=0; j<4; j++)
      EXPECT_NEAR( R(i,j), trueR(i,j), 1e-6);

  int shift=rank_*localSize_;
  for (auto i=0; i<localSize_; i++)
    for (auto j=0; j<Q.localNumVectors(); j++)
      EXPECT_NEAR( Q(i,j), trueQ(i+shift,j), 1e-6);

}
