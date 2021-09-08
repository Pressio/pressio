#include "pressio/ops.hpp"
#include "pressio/rom/tpw_roms/tpql_operators.hpp"
#include "pressio/rom/tpw_roms/tpql_rom_datatypes.hpp"

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

  using execution_space = Kokkos::DefaultExecutionSpace;
  using kll   = Kokkos::LayoutLeft;
  mutable Kokkos::View<scalar_t*, kll, execution_space> params_; 

  MyTpetraApp(rcpcomm_t comm): comm_{comm}, params_("params",3){
    params_(0) = 2;
    params_(1) = 7;
    params_(2) = 4;
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
    f_v[0] = y_v[0]*y_v[1]*y_v[1] + params_(0)*params_(1)*params_(1);
    f_v[1] = y_v[1]*y_v[0] + y_v[0]*y_v[0] + params_(2)*params_(2)*y_v[1]*params_(1);
  }
  
  state_t const & initialCondition() const{
    return *U0_;
  }

  template <typename params_t>
  void updateScalarParameters(params_t & params_in) const {
    for (int i = 0; i < 3; i++){
      params_(i) = params_in(i);
    }
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
  mutable Eigen::Matrix<double,-1,1> params_;

  MyEigenApp() : params_(3){
    params_(0) = 2;
    params_(1) = 7;
    params_(2) = 4;
  }

  template <typename params_t>
  void updateScalarParameters(params_t & params_in) const{
    for (int i = 0; i < 3; i++){
      params_(i) = params_in(i);
    }
  }

  
  Eigen::Matrix<double,-1,1> createVelocity() const{
    velocity_t f(N_);
    return f;
  }
  
  template <typename state_t, typename velocity_t> 
  void velocity(const state_t & y, const double t, velocity_t & f)
  const {
    f(0) = y(0)*y(1)*y(1) + params_(0)*params_(1)*params_(1) ;
    f(1) = y(1)*y(0) + y(0)*y(0) + params_(2)*params_(2)*y(1)*params_(1);
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
  using eigen_data_t = typename ::pressio::rom::experimental::RomDataTypeEigen<scalar_t>;
  using kokkos_data_t = typename ::pressio::rom::experimental::RomDataTypeKokkos<scalar_t>;
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
 
    Eigen::Matrix<scalar_t , -1,1> mu(3);
    mu(0) = 2;
    mu(1) = 7;
    mu(2) = 4;

    //auto romDataTypeEnum = pressio::rom::experimental::RomDataType::Kokkos;
    auto PhiTJPhiKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasisKokkos<kokkos_data_t>(appObj, u, mu, 0., Phi);
    auto PhiTJPhiEigen = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasisEigen<eigen_data_t>(appObj, u, mu, 0., Phi);

    double tol = 1e-3; 
    Eigen::Matrix<scalar_t, -1,-1> PhiTJPhi_exact(N,N);
    PhiTJPhi_exact(0,0) = 4;
    PhiTJPhi_exact(0,1) = 12;
    PhiTJPhi_exact(1,0) = 8;
    PhiTJPhi_exact(1,1) = 115;

    for (int i=0; i < 2; i++){
      for (int j=0; j < 2; j++){
        auto e1 = PhiTJPhiKokkos(i,j) - PhiTJPhi_exact(i,j);
        auto e2 = PhiTJPhiEigen(i,j) - PhiTJPhi_exact(i,j);
        if (std::abs(e1) > tol){
          std::cout << " error in Jacobian routine for Kokkos, Kokkos evaluates to " <<  PhiTJPhiKokkos(i,j)  << " vs " << PhiTJPhi_exact(i,j) << std::endl;}
        if (std::abs(e2) > tol){
          std::cout << " error in Jacobian routine for Eigen, Eigen evaluates to " <<  PhiTJPhiEigen(i,j)  << " vs " << PhiTJPhi_exact(i,j) << std::endl;}
      }
    }
    //std::cout << PhiTJPhi - PhiTJPhi_exact << std::endl;

    auto PhiTHJPhiKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasisKokkos<kokkos_data_t>(appObj, u, mu, 0., Phi);
    auto PhiTHJPhiEigen = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasisEigen<eigen_data_t>(appObj, u, mu, 0., Phi);

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
          auto e1 =  PhiTHJPhiKokkos(i,j,k) - Hexact_list(indx);
          auto e2 = PhiTHJPhiEigen[i][j][k] - Hexact_list(indx);
          if (std::abs(e1) > tol){
            std::cout << " error in Hessian routine for Kokkos, Kokkos evaluates to " <<  PhiTHJPhiKokkos(i,j,k)  << " vs " << Hexact_list(indx) << std::endl;}
          if (std::abs(e2) > tol){
            std::cout << " error in Hessian routine for Eigen, Eigen evaluates to " <<  PhiTHJPhiEigen[i][j][k]  << " vs " << Hexact_list(indx) << std::endl;}
          indx += 1;}}}


    auto PhiTJParamsKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobianKokkos<kokkos_data_t>(appObj, u, mu, 0., Phi);
    auto PhiTJParamsEigen = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobianEigen<eigen_data_t>(appObj, u, mu, 0., Phi);

    Eigen::Matrix<scalar_t, -1,-1> PhiTJParams_exact(N,3);
    PhiTJParams_exact(0,0) = 49;
    PhiTJParams_exact(0,1) = 28;
    PhiTJParams_exact(0,2) = 0;

    PhiTJParams_exact(1,0) = 0;
    PhiTJParams_exact(1,1) = 32;
    PhiTJParams_exact(1,2) = 112;

    for (int i=0; i < 2; i++){
      for (int j=0; j < 3; j++){
        auto e1 = PhiTJParamsKokkos(i,j) - PhiTJParams_exact(i,j);
        auto e2 = PhiTJParamsEigen(i,j) - PhiTJParams_exact(i,j);
        if (std::abs(e1) > tol){
          std::cout << " error in parameter Jac routine for Kokkos, Kokkos evaluates to " <<  PhiTJParamsKokkos(i,j)  << " vs " << PhiTJParams_exact(i,j) << std::endl;}
        if (std::abs(e2) > tol){
          std::cout << " error in parameter Jac routine for Eigen, Eigen evaluates to " <<  PhiTJParamsEigen(i,j)  << " vs " << PhiTJParams_exact(i,j) << std::endl;}
      }
    }


    

    auto PhiTHParamsKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessianKokkos<kokkos_data_t>(appObj, u, mu, 0., Phi);
    auto PhiTHParamsEigen = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessianEigen<eigen_data_t>(appObj, u, mu, 0., Phi);

    kokkos_data_t::dense_hessian_t PhiTHParams_exact("Hparams",2,3,3);
    PhiTHParams_exact(0,0,0) = 0;
    PhiTHParams_exact(0,0,1) = 14.;
    PhiTHParams_exact(0,0,2) = 0;
    PhiTHParams_exact(0,1,0) = 14.;
    PhiTHParams_exact(0,1,1) = 4;
    PhiTHParams_exact(0,1,2) = 0;
    PhiTHParams_exact(0,2,0) = 0;
    PhiTHParams_exact(0,2,1) = 0;
    PhiTHParams_exact(0,2,2) = 0;

    PhiTHParams_exact(1,0,0) = 0;
    PhiTHParams_exact(1,0,1) = 0;
    PhiTHParams_exact(1,0,2) = 0;
    PhiTHParams_exact(1,1,0) = 0;
    PhiTHParams_exact(1,1,1) = 0;
    PhiTHParams_exact(1,1,2) = 16.;
    PhiTHParams_exact(1,2,0) = 0;
    PhiTHParams_exact(1,2,1) = 16.;
    PhiTHParams_exact(1,2,2) = 28.;

    for (int i=0; i < 2; i++){
      for (int j=0; j < 3; j++){
        for (int k=0; k < 3; k++){
          auto e1 = PhiTHParamsKokkos(i,j,k) - PhiTHParams_exact(i,j,k);
          auto e2 = PhiTHParamsEigen[i][j][k] - PhiTHParams_exact(i,j,k);
          if (std::abs(e1) > tol){
            std::cout << " error in parameter Hessian routine for Kokkos, Kokkos evaluates to " <<   PhiTHParamsKokkos(i,j,k)  << " vs " << PhiTHParams_exact(i,j,k) << std::endl;}
          if (std::abs(e2) > tol){
            std::cout << " error in parameter Hessian routine for Eigen, Eigen evaluates to " <<  PhiTHParamsEigen[i][j][k]  << " vs " << PhiTHParams_exact(i,j,k)  << std::endl;}
        }
      }
    }


    auto PhiTHMixedKokkos = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessianKokkos<kokkos_data_t>(appObj, u, mu, 0., Phi);
    auto PhiTHMixedEigen = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessianEigen<eigen_data_t>(appObj, u, mu, 0., Phi);

    kokkos_data_t::dense_hessian_t PhiTHMixed_exact("Hparams",2,2,3);
    PhiTHMixed_exact(1,1,1) = 16.;
    PhiTHMixed_exact(1,1,2) = 56.;

    for (int i=0; i < 2; i++){
      for (int j=0; j < 2; j++){
        for (int k=0; k < 3; k++){
          auto e1 = PhiTHMixedKokkos(i,j,k) - PhiTHMixed_exact(i,j,k);
          auto e2 = PhiTHMixedEigen[i][j][k] - PhiTHMixed_exact(i,j,k);
          if (std::abs(e1) > tol){
            std::cout << " error in parameter mixed Hessian routine for Kokkos, Kokkos evaluates to " <<   PhiTHMixedKokkos(i,j,k)  << " vs " << PhiTHMixed_exact(i,j,k) << std::endl;}
          if (std::abs(e2) > tol){
            std::cout << " error in parameter mixed Hessian routine for Eigen, Eigen evaluates to " <<  PhiTHMixedEigen[i][j][k]  << " vs " << PhiTHMixed_exact(i,j,k)  << std::endl;}
        }
      }
    }


  /*

  */
  }

  return 0;
}
