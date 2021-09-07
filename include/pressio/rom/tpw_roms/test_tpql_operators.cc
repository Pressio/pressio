#include "pressio/ops.hpp"
#include "pressio/rom/tpw_roms/tpql_operators.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>



class MyTpetraApp
{
protected:
  using map_t		= Tpetra::Map<>;
  using nativeVec	= Tpetra::Vector<>;

  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  int N_ = 2;
  int myRank_{};
  int totRanks_{};
  template<typename T> using stdrcp = std::shared_ptr<T>;

public:
  using velocity_t = Tpetra::Vector<>;
  using state_t = Tpetra::Vector<>;
  using scalar_t = double;

  MyTpetraApp(rcpcomm_t comm): comm_{comm}{
    this->setup();
  }

  rcpmap_t getDataMap(){
    return dataMap_;
  };

  
  Tpetra::Vector<> createVelocity() const{
    velocity_t f(dataMap_);
    return f;
  }
  
  template <typename state_t, typename velocity_t> 
  void velocity(const state_t & y, const double t, velocity_t & f)
  const {
    auto f_v = f.getDataNonConst();
    auto y_v = y.getData();
    f_v[0] = y_v[0]*y_v[1]*y_v[1];
    f_v[1] = y_v[1]*y_v[0] + y_v[0]*y_v[0];
  }
  
  state_t const & initialCondition() const{
    return *U0_;
  }

protected:
  void setup()
  {
    //myRank_ =  comm_->getRank();
    //totRanks_ =  comm_->getSize();
    // distribute cells
    dataMap_ = Teuchos::rcp(new map_t(2, 0, comm_));
    U0_ = std::make_shared<state_t>(dataMap_);
    auto ud = (*U0_).getDataNonConst();
    ud[0] = 3.;
    ud[1] = 2.;
//    U0_->putScalar(1.0);

  }
  Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
  rcpcomm_t comm_{};
  rcpmap_t dataMap_{};
  mutable stdrcp<state_t> U0_{}; // initial state vector

};


class MyEigenApp
{
public:
  MyEigenApp(){}

  
  Eigen::Matrix<double,-1,1> createVelocity() const{
    velocity_t f(N_);
    return f;
  }
  
  template <typename state_t, typename velocity_t> 
  void velocity(const state_t & y, const double t, velocity_t & f)
  const {
    f(0) = y(0)*y(1)*y(1);
    f(1) = y(1)*y(0) + y(0)*y(0);
  }
  
  Eigen::Matrix<double,-1,1> initialCondition() const{
    state_t y(N_);
    y(0) = 3.;
    y(1) = 2.;
    return y;
  }


private: 
  int N_ = 2;
  using velocity_t = Eigen::Matrix<double,-1,1>;
  using state_t = Eigen::Matrix<double,-1,1>;

};


int main(int argc, char *argv[])
{
  using scalar_t = double;



  int N = 2;


  // Eigen tests
  /*
  auto appObj = MyEigenApp();
  Eigen::Matrix<scalar_t, -1,-1> Phi(N,N);
  Phi.setZero();
  Phi(0,0) = 1;
  Phi(1,1) = 1;
  auto u = appObj.initialCondition();
 
  Eigen::Matrix<scalar_t , -1,1> mu(1);
  auto PhiTJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasis(appObj, u, mu, 0., Phi);

  Eigen::Matrix<scalar_t, -1,-1> PhiTJPhi_exact(N,N);
  PhiTJPhi_exact(0,0) = 4;
  PhiTJPhi_exact(0,1) = 12;
  PhiTJPhi_exact(1,0) = 8;
  PhiTJPhi_exact(1,1) = 3;


  auto PhiTHJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasis(appObj, u, mu, 0., Phi);

  Eigen::Matrix<scalar_t, -1, 1> Hexact_list(8);
  Hexact_list(0) = 0;
  Hexact_list(1) = 4;
  Hexact_list(2) = 4;
  Hexact_list(3) = 6;
  Hexact_list(4) = 2;
  Hexact_list(5) = 1;
  Hexact_list(6) = 1;
  Hexact_list(7) = 0;

  int indx = 0;
  for (int i=0; i < 2; i++){
    for (int j=0; j < 2; j++){
      for (int k=0; k < 2 ; k++){
        std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTHJPhi[i][j][k] - Hexact_list(indx) << std::endl;
        indx += 1;}}}


  auto PhiTParamJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobian(appObj, u, mu, 0., Phi);
  std::cout << PhiTParamJPhi << std::endl;

  auto PhiTParamHPhi = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessian(appObj, u, mu, 0., Phi);
  indx = 0;
  for (int i=0; i < 2; i++){
    for (int j=0; j < 1; j++){
      for (int k=0; k < 1 ; k++){
        std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTParamHPhi[i][j][k]  << std::endl;
        indx += 1;}}}

  */

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  { 
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rcpcomm_t Comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    auto appObj = MyTpetraApp(Comm);
    auto dataMap = appObj.getDataMap();
    using basis_t	= Tpetra::MultiVector<>;
    basis_t Phi(dataMap,N,true);
    auto Phi0 = Phi.getDataNonConst(0);
    Phi0[0] = 1.;
    auto Phi1 = Phi.getDataNonConst(1);
    Phi1[1] = 1.;

    auto u = appObj.initialCondition();
 
    Eigen::Matrix<scalar_t , -1,1> mu(1);
    auto romDataTypeEnum = pressio::rom::experimental::RomDataType::Kokkos;
    auto PhiTJPhiKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasisKokkos(appObj, u, mu, 0., Phi);
    auto PhiTJPhiEigen = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasisEigen(appObj, u, mu, 0., Phi);

    
    Eigen::Matrix<scalar_t, -1,-1> PhiTJPhi_exact(N,N);
    PhiTJPhi_exact(0,0) = 4;
    PhiTJPhi_exact(0,1) = 12;
    PhiTJPhi_exact(1,0) = 8;
    PhiTJPhi_exact(1,1) = 3;

    for (int i=0; i < 2; i++){
      for (int j=0; j < 2; j++){
        std::cout << PhiTJPhiKokkos(i,j) <<  std::endl;}}
    std::cout << PhiTJPhiEigen <<  std::endl;
    std::cout << PhiTJPhi_exact <<  std::endl;

    //std::cout << PhiTJPhi - PhiTJPhi_exact << std::endl;

    auto PhiTHJPhiKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasisKokkos(appObj, u, mu, 0., Phi);
    auto PhiTHJPhiEigen = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasisEigen(appObj, u, mu, 0., Phi);

    Eigen::Matrix<scalar_t, -1, 1> Hexact_list(8);
    Hexact_list(0) = 0;
    Hexact_list(1) = 4;
    Hexact_list(2) = 4;
    Hexact_list(3) = 6;
    Hexact_list(4) = 2;
    Hexact_list(5) = 1;
    Hexact_list(6) = 1;
    Hexact_list(7) = 0;

    int indx = 0;
    for (int i=0; i < 2; i++){
      for (int j=0; j < 2; j++){
        for (int k=0; k < 2 ; k++){
//          std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTHJPhi[i][j][k] - Hexact_list(indx) << std::endl;
          std::cout << "i = " << i << " j = " << j << " k = " << k << " sol = " << PhiTHJPhiKokkos(i,j,k) - Hexact_list(indx) << std::endl;
          indx += 1;}}}


    auto PhiTJParamsKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobianKokkos(appObj, u, mu, 0., Phi);
    auto PhiTJParamsEigen = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobianEigen(appObj, u, mu, 0., Phi);

    auto PhiTHParamsKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessianKokkos(appObj, u, mu, 0., Phi);
    auto PhiTHParamsEigen = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessianEigen(appObj, u, mu, 0., Phi);

    auto PhiTHMixedKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessianKokkos(appObj, u, mu, 0., Phi);
    auto PhiTHMixedEigen = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessianEigen(appObj, u, mu, 0., Phi);


    /*
  */
  }

  return 0;
}
