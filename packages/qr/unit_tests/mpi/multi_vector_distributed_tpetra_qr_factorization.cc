
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "../../src/qr_tpetra_multi_vector.hpp"

struct tpetraMultiVectorR9C4Fixture
  : public ::testing::Test{

public:
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using mvec_t = Tpetra::MultiVector<>;
  using ST = typename mvec_t::scalar_type;
  using LO = typename mvec_t::local_ordinal_type;
  using GO = typename mvec_t::global_ordinal_type;

  int rank_;
  int numProc_;
  const int localSize_ = 3;
  const int numVecs_ = 4;
  int numGlobalEntries_;
  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> contigMap_;
  mvec_t * mv_;

  virtual void SetUp(){
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    comm_ = Teuchos::rcp (new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    rank_ = comm_->getRank();
    numProc_ = comm_->getSize();
    assert(numProc_==3);
    numGlobalEntries_ = numProc_ * localSize_;
    contigMap_ = Teuchos::rcp(new map_t(numGlobalEntries_,0,comm_));
    mv_ = new mvec_t(contigMap_, numVecs_);
  }

  virtual void TearDown(){
    delete mv_;
  }
};
//-----------------------------------------------------------


TEST_F(tpetraMultiVectorR9C4Fixture,
       TpetraMultiVectorQRFactorization){
  using namespace rompp;

  // --------------------------------------------
  // construct and fill multivector wrapper
  // --------------------------------------------
  using nat_mv_t = typename tpetraMultiVectorR9C4Fixture::mvec_t;
  using mvec_t = core::MultiVector<nat_mv_t>;
  using mv_device_t = typename core::details::traits<mvec_t>::device_t;

  mvec_t MV( *mv_ );
  MV.setZero();
  // get trilinos tpetra multivector object
  auto trilD = MV.data();
  trilD->sync<Kokkos::HostSpace>();

  /*--------------------------------------------
   * (1): modify the host view and then sync
   * most likely, host and device will be same unless we run CUDA
   * so in theory we should not worry about syncing but it
   * does not hurt to do it anyway
  //--------------------------------------------*/
  auto v2d = trilD->getLocalView<Kokkos::HostSpace>();
  auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
  //we are going to change the host view
  trilD->modify<Kokkos::HostSpace>();

  if(rank_==0){
    v2d(0,0) = 3.2; v2d(0,1) = 1.2;  v2d(0,2) = 1.;
    v2d(1,0) = 1.2; v2d(1,2) = -2.2;
    v2d(2,1) = 4.0; v2d(2,3) = -2.;
  }
  if(rank_==1){
    v2d(0,1) = 4.;
    v2d(1,2) = -1.; v2d(1,3) = -4.;
    v2d(2,0) = 0.2; v2d(2,1) = 5.;  v2d(2,2) = 1.;
  }
  if(rank_==2){
    v2d(0,0) = 1.; v2d(0,1) = 1.1; v2d(0,2) = 1.25; v2d(0,3) = -3.;
    v2d(1,2) = 1.; v2d(1,1) = 0.1111; v2d(1,3) = 6.;
  }
  // sync from host to device
  trilD->sync<mv_device_t>();

  // do QR
  using R_type = rompp::core::Matrix<Eigen::MatrixXd>;
  rompp::qr::QRSolver<mvec_t, rompp::core::MultiVector, R_type> qrObj;
  qrObj.compute(MV);
  const auto & Q = qrObj.cRefQFactor();
  const auto & R = qrObj.cRefRFactor();
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
  for (auto j=0; j<Q.localNumVectors(); j++){
    auto colData = Q.data()->getData(j);
    for (auto i=0; i<localSize_; i++){
      EXPECT_NEAR( colData[i], trueQ(i+shift,j), 1e-6);
    }
  }

}
