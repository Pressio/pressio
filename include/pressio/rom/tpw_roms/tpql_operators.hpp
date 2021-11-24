#include "pressio/ops.hpp"
#include "tpql_rom_datatypes.hpp"
#ifndef TPQL_OPS_HPP_
#define TPQL_OPS_HPP_

namespace pressio{ namespace rom{ namespace experimental{



template <typename scalar_t>
auto getCol(const Tpetra::BlockMultiVector<> & Phi,const int i){

  using map_t   = Tpetra::Map<>;

  using go_t    = typename map_t::global_ordinal_type;
  using lo_t    = typename map_t::local_ordinal_type;

  auto meshMap = Phi.getPointMap();
  auto PhiMultiVector = Phi.getMultiVectorView(); 
  auto PhiCol = PhiMultiVector.getVector(i);
  Tpetra::BlockVector<> PhiColBlockVector(*PhiCol , meshMap, 1);
  return PhiColBlockVector;
}


template <typename scalar_t>
auto getCol(const Tpetra::MultiVector<> & Phi,const int i){
  auto PhiCol = Phi.getVector(i);
  return *PhiCol;
}



template <typename scalar_t, typename matrix_t>
Eigen::Matrix<scalar_t, -1 , 1> getCol(const matrix_t & Phi,const int i){
  Eigen::Matrix<scalar_t,-1,1> PhiCol = Phi.col(i);
  return PhiCol;
}

template <typename scalar_t>
void setCol(Eigen::Matrix<scalar_t,-1,-1> & Phi,const Eigen::Matrix<scalar_t,-1,1>  & PhiCol,const int i){
  Phi.col(i) = PhiCol;
}


template <typename scalar_t, typename layout_t, typename exec_space_t>
void setCol(Kokkos::View<scalar_t**,layout_t,exec_space_t> & Phi,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & PhiCol, const int i){
 auto numRows = ::pressio::ops::extent(Phi,0);
 for (int j = 0; j < numRows; j++){
   Phi(j,i) = PhiCol(j);
 }
}


template <typename scalar_t>
void setHessianCol(std::vector<std::vector<std::vector<scalar_t> > > & H , const Eigen::Matrix<scalar_t,-1,1>  & HCol,const int i, const int j )
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H[k][i][j] = HCol(k);
  } 
}

template <typename scalar_t, typename layout_t, typename exec_space_t>
void setHessianCol(Kokkos::View<scalar_t***,layout_t,exec_space_t> & H,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & HCol,const int i,const int j)
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H(k,i,j) = HCol(k);
  } 
}

template <typename scalar_t>
void setSymmetricHessianCol(std::vector<std::vector<std::vector<scalar_t> > > & H , const Eigen::Matrix<scalar_t,-1,1>  & HCol,const int i, const int j )
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H[k][i][j] = HCol(k);
    H[k][j][i] = HCol(k);
  } 
}

template <typename scalar_t, typename layout_t, typename exec_space_t>
void setSymmetricHessianCol(Kokkos::View<scalar_t***,layout_t,exec_space_t> & H,const Kokkos::View<scalar_t*,layout_t,exec_space_t> & HCol,const int i,const int j)
{
  auto dim = ::pressio::ops::extent(HCol,0);
  for (int k = 0; k < dim; k++){
    H(k,i,j) = HCol(k);
    H(k,j,i) = HCol(k);
  } 
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename dense_matrix_t, typename vector_t>
void computeTestBasisTransposeTimesMixedHessian_calc(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi,const basis_t & TestBasis, dense_matrix_t & BasisTransposeTimesMixedHessian, vector_t & BasisTransposeTimesMixedHessianCol, const scalar_t epsilon1, const scalar_t epsilon2)

{
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto y0_norm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  auto yp = ::pressio::ops::clone(y);
  auto muP = ::pressio::ops::clone(mu);

  //scalar_t epsState = 1e-3;
  scalar_t epsParams = epsilon2;

  //scalar_t scale = 0.25/(epsState*epsParams);

  for  (int i = 0; i < romDim; i++){
    auto PhiCol = getCol<scalar_t>(Phi,i);
    auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
    scalar_t epsState = std::sqrt( (1. + y0_norm)*epsilon1)/PhiColNorm;
    scalar_t scale = 0.25 / (epsState * epsParams);
    for (int j=0; j < numParams; j++){
      ::pressio::ops::set_zero(workingVector);
      ::pressio::ops::update(yp,0.,y,1.,PhiCol,epsState);
      muP(j) += epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);
      muP(j) -= epsParams;

      muP(j) -= epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
      muP(j) += epsParams;

      ::pressio::ops::update(yp,1.,PhiCol,-2.*epsState);
      muP(j) += epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
      muP(j) -= epsParams;

      muP(j) -= epsParams;
      appObj.updateScalarParameters(muP);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);
      muP(j) += epsParams;

      ::pressio::ops::product(::pressio::transpose(),1., TestBasis ,workingVector , 0.,BasisTransposeTimesMixedHessianCol);
      setHessianCol(BasisTransposeTimesMixedHessian,BasisTransposeTimesMixedHessianCol,i,j);
    }
  }
  appObj.updateScalarParameters(mu);
}




template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesMixedHessianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const basis_t & TestBasis, const scalar_t epsilon1, const scalar_t epsilon2)

{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesMixedHessian("PhiTH_params", romDim,romDim,numParams);
  vector_t HCol("HCol",romDim);
  computeTestBasisTransposeTimesMixedHessian_calc(appObj,y,mu,t,Phi,TestBasis,BasisTransposeTimesMixedHessian,HCol,epsilon1, epsilon2);
  return BasisTransposeTimesMixedHessian; 
}





template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
std::vector< std::vector< std::vector< scalar_t > > > computeTestBasisTransposeTimesMixedHessianEigen(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const basis_t & TestBasis, const scalar_t epsilon1, const scalar_t epsilon2)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t = typename rom_data_t::vector_t;
  using dense_hessian_t = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesMixedHessian(romDim);
  vector_t BasisTransposeTimesMixedHessianCol(romDim);
  BasisTransposeTimesMixedHessianCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesMixedHessian[i].resize(romDim);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < romDim; j++){
      BasisTransposeTimesMixedHessian[i][j].resize(numParams);
      for (int k=0; k < numParams; k++){
        BasisTransposeTimesMixedHessian[i][j][k] = 0.;
      }
    }
  }

  computeTestBasisTransposeTimesMixedHessian_calc(appObj,y,mu,t,Phi,TestBasis,BasisTransposeTimesMixedHessian,BasisTransposeTimesMixedHessianCol,epsilon1,epsilon2);
  return BasisTransposeTimesMixedHessian; 
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeTestBasisTransposeTimesParameterHessian_calc(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi,const basis_t & TestBasis,  hessian_t & BasisTransposeTimesParameterHessian, vector_t & BasisTransposeTimesParameterHessianCol, const scalar_t epsilon)
{
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto numParams = ::pressio::ops::extent(mu,0);
  scalar_t eps = epsilon;

  auto muP = ::pressio::ops::clone(mu);

  scalar_t scale1 = 1./(12.*eps*eps);
  scalar_t scale2 = 0.25/(eps*eps);

  for  (int i = 0; i < numParams; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        muP(i) += 2.*eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.*scale1);
        muP(i) -= 2.*eps;

        muP(i) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale1);
        muP(i) -= eps;

        ::pressio::ops::update(workingVector,1.,f0,-30.*scale1);
         
        muP(i) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale1);
        muP(i) += eps;

        muP(i) -= 2.*eps; 
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale1);
        muP(i) += 2.*eps;
      }
      else{
        muP(i) += eps;
        muP(j) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.*scale2);
        muP(i) -= eps;
        muP(j) -= eps;

        muP(i) += eps;
        muP(j) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale2);
        muP(i) -= eps;
        muP(j) += eps;

        muP(i) -= eps;
        muP(j) += eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale2);
        muP(i) += eps;
        muP(j) -= eps;

        muP(i) -= eps;
        muP(j) -= eps;
        appObj.updateScalarParameters(muP);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(y,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.*scale2);
        muP(i) += eps;
        muP(j) += eps;
      }
      ::pressio::ops::product(::pressio::transpose(),1., TestBasis ,workingVector , 0.,BasisTransposeTimesParameterHessianCol);


      setSymmetricHessianCol(BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,i,j);
    }
  }
  //reset parameters
  appObj.updateScalarParameters(mu);
}



template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesParameterHessianEigen(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const basis_t & TestBasis, const scalar_t epsilon)
{

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);
  using vector_t = typename rom_data_t::vector_t;
  using dense_hessian_t = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesParameterHessian(romDim);
  vector_t BasisTransposeTimesParameterHessianCol(romDim);
  BasisTransposeTimesParameterHessianCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesParameterHessian[i].resize(numParams);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < numParams; j++){
      BasisTransposeTimesParameterHessian[i][j].resize(numParams);
      for (int k=0; k < numParams; k++){
        BasisTransposeTimesParameterHessian[i][j][k] = 0.;
      }
    }
  }
  computeTestBasisTransposeTimesParameterHessian_calc(appObj,y,mu,t,Phi,TestBasis, BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,epsilon);
  return BasisTransposeTimesParameterHessian; 
}

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesParameterHessianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const  scalar_t t,const basis_t & Phi, const basis_t & TestBasis, const scalar_t epsilon)
{

  auto romDim = ::pressio::ops::extent(Phi,1);
  auto numParams = ::pressio::ops::extent(mu,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    = typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesParameterHessian("PhiTH_params", romDim,numParams,numParams);
  vector_t BasisTransposeTimesParameterHessianCol("HCol",romDim);
  computeTestBasisTransposeTimesParameterHessian_calc(appObj,y,mu,t,Phi,TestBasis,BasisTransposeTimesParameterHessian,BasisTransposeTimesParameterHessianCol,epsilon);
  return BasisTransposeTimesParameterHessian; 
}


/*--------------------------------------------------
Phi^T J_params computation
--------------------------------------------------*/
template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename dense_matrix_t, typename vector_t>
void computeTestBasisTransposeTimesParameterJacobian_calc(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const basis_t & TestBasis, dense_matrix_t & PhiTJPhi, vector_t & PhiTJPhiCol, const scalar_t epsilon)
{
  auto ftmp = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  auto workingVector = appObj.createVelocity();

  auto numParams = ::pressio::ops::extent(mu,0);
  auto fomDim = ::pressio::ops::extent(y,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  scalar_t eps = epsilon;
  scalar_t scale = 0.5/eps;
  auto muP = ::pressio::ops::clone(mu);
  for  (int i = 0; i < numParams; i++){
    ::pressio::ops::set_zero(workingVector);

    muP(i) += eps;
    appObj.updateScalarParameters(muP); 
    ::pressio::ops::set_zero(ftmp);
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);

    muP(i) -= 2.*eps;
    ::pressio::ops::set_zero(ftmp);
    appObj.updateScalarParameters(muP); 
    appObj.velocity(y,t,ftmp);
    ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

    ::pressio::ops::product(::pressio::transpose(),1.,TestBasis,workingVector , 0., PhiTJPhiCol);
    setCol(PhiTJPhi,PhiTJPhiCol,i);
  
  }
  // reset parameters to original state
  appObj.updateScalarParameters(mu);
} 


template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesParameterJacobianEigen(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const basis_t & TestBasis, const scalar_t epsilon)
{
  auto numParams = ::pressio::ops::extent(mu,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  using dense_matrix_t = typename rom_data_t::dense_matrix_t;
  using vector_t = typename rom_data_t::vector_t;

  dense_matrix_t PhiTJPhi(romDim,numParams);
  PhiTJPhi.setZero();

  vector_t PhiTJPhiCol(romDim);
  PhiTJPhiCol.setZero();

  computeTestBasisTransposeTimesParameterJacobian_calc(appObj,y,mu,t,Phi,TestBasis, PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesParameterJacobianKokkos(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,const basis_t & TestBasis, const scalar_t epsilon)
{
  auto numParams = ::pressio::ops::extent(mu,0);
  auto romDim = ::pressio::ops::extent(Phi,1);

  using dense_matrix_t = typename rom_data_t::dense_matrix_t;
  using vector_t = typename rom_data_t::vector_t;

  dense_matrix_t PhiTJPhi("PhiJPphi_params",romDim,numParams);
  vector_t PhiTJPhiCol("PhiJPphiCol_params",romDim);

  computeTestBasisTransposeTimesParameterJacobian_calc(appObj,y,mu,t,Phi,TestBasis,PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 

//-----------------------------------



/*---------------------------------
Phi^T f computation
----------------------------------*/
template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
Eigen::Matrix<scalar_t,-1,1> computeTestBasisTransposeTimesVelocityEigen(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & TestBasis)
{
  auto romDim = ::pressio::ops::extent(TestBasis,1);

  // Create temporary working vectors
  auto f = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity(y,t,f);

  using vector_t = typename rom_data_t::vector_t;
  vector_t PhiTfCol(romDim);
  PhiTfCol.setZero();

  ::pressio::ops::product(::pressio::transpose(),1., TestBasis,f, 0.,PhiTfCol);

  auto JPhiVecViewA = TestBasis.getMultiVectorView();
  //std::cout << " hereB " <<  JPhiVecViewA.getLocalLength() << std::endl;

  return PhiTfCol;
} 

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesVelocityKokkos(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & TestBasis)
{
  auto romDim = ::pressio::ops::extent(TestBasis,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity(y,t,f);

  using vector_t		= typename rom_data_t::vector_t;
  vector_t PhiTfCol("PhiTfCol",romDim);

  ::pressio::ops::product(::pressio::transpose(),1., TestBasis,f, 0.,PhiTfCol);
  return PhiTfCol;
} 


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename jphi_t>
void computeJacobianTimesBasis_calc(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,jphi_t & JPhi, const scalar_t epsilon){
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto ynorm = ::pressio::ops::norm2(y);
  // Create temporary working vectors
  auto fp = appObj.createVelocity();
  auto yp = ::pressio::ops::clone(y);
  auto workingVector = appObj.createVelocity();
  appObj.updateScalarParameters(mu);

  auto JPhi_mv = JPhi.getMultiVectorView();

  for  (int i = 0; i < romDim; i++){
    ::pressio::ops::set_zero(workingVector);

    auto PhiCol = getCol<scalar_t>(Phi,i);
    auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
    scalar_t eps = std::sqrt( (1. + ynorm)*epsilon ) / PhiColNorm;
    scalar_t scale = 0.5/eps;

    ::pressio::ops::update(yp,0.,y,1.,PhiCol,eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,0.,fp,1.*scale);

    ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,1.,fp,-1.*scale);

    auto mode_view = JPhi_mv.getVectorNonConst(i);
    auto workingVector_view = workingVector.getVectorView();
    mode_view->update(1., workingVector_view, 0.0);
  }
}


template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
void computeJacobianTimesBasis(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi,const basis_t & JPhi, const scalar_t epsilon)
{
  //auto JPhi = appObj.createApplyJacobianResult(Phi);
  computeJacobianTimesBasis_calc(appObj,y,mu,t,Phi,JPhi,epsilon);
} 


// 
template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename test_basis_t, typename dense_matrix_t, typename vector_t>
void computeTestBasisTransposeTimesJacobianTimesBasis_calc(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const test_basis_t & TestBasis, dense_matrix_t & PhiTJPhi,vector_t & PhiTJPhiCol, const scalar_t epsilon){

  auto JPhiVecViewA = TestBasis.getMultiVectorView();
  //std::cout << " JPhiLocalLength in computeTestBasisTransposeTimesJacobianTimesBasis_calc " <<  JPhiVecViewA.getLocalLength() << std::endl;


  auto romDim = ::pressio::ops::extent(PhiTJPhi,0);

  auto ynorm = ::pressio::ops::norm2(y);
  // Create temporary working vectors
  auto fp = appObj.createVelocity();
  auto yp = ::pressio::ops::clone(y);
  auto workingVector = appObj.createVelocity();
  appObj.updateScalarParameters(mu);

  for  (int i = 0; i < romDim; i++){
    ::pressio::ops::set_zero(workingVector);

    auto PhiCol = getCol<scalar_t>(Phi,i);
    auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
    scalar_t eps = std::sqrt( (1. + ynorm)*epsilon ) / PhiColNorm;
    scalar_t scale = 0.5/eps;
    //std::cout << " a " << std::endl;

    ::pressio::ops::update(yp,0.,y,1.,PhiCol,eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,0.,fp,1.*scale);
    //std::cout << " b " << std::endl;

    ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
    appObj.velocity(yp,t,fp);
    ::pressio::ops::update(workingVector,1.,fp,-1.*scale);
    //std::cout << " c " << std::endl;

    //std::cout << "JPhi Size " << ::pressio::ops::extent(TestBasis,0) << " " << ::pressio::ops::extent(TestBasis,1) << std::endl;
    //std::cout << "workinv vector size " << ::pressio::ops::extent(workingVector,0) << std::endl;

    auto JPhiVecView = TestBasis.getMultiVectorView();
    auto PhiVecView = Phi.getMultiVectorView();

    auto workingVectorView = workingVector.getMultiVectorView();
    //std::cout << " JPhiLocalLength in computeTestBasisTransposeTimesJacobianTimesBasis_calc loop " <<  JPhiVecView.getLocalLength() << std::endl;
    //std::cout << " PhiLocalLength in computeTestBasisTransposeTimesJacobianTimesBasis_calc loop " <<  PhiVecView.getLocalLength() << std::endl;
    //std::cout << " working vector local Length in computeTestBasisTransposeTimesJacobianTimesBasis_calc loop " <<  workingVectorView.getLocalLength() << std::endl;

    ::pressio::ops::product(::pressio::transpose(),1., TestBasis,workingVector , 0., PhiTJPhiCol);

    //std::cout << " d " << std::endl;

    setCol(PhiTJPhi,PhiTJPhiCol,i);


  }
}

template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename test_basis_t>
auto computeTestBasisTransposeTimesJacobianTimesBasisKokkos(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const  basis_t & Phi, const test_basis_t & TestBasis, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);

  using vector_t		= typename rom_data_t::vector_t;
  using dense_matrix_t		= typename rom_data_t::dense_matrix_t;

  dense_matrix_t PhiTJPhi("PhiTJPhi",romDim,romDim); 
  vector_t PhiTJPhiCol("PhiTJPhiCol",romDim);
  computeTestBasisTransposeTimesJacobianTimesBasis_calc(appObj,y,mu,t,Phi,TestBasis,PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 


template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename test_basis_t>
auto computeTestBasisTransposeTimesJacobianTimesBasisEigen(const app_t & appObj, state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,const test_basis_t & TestBasis, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);

  // Create result
  using vector_t		= typename rom_data_t::vector_t;
  using dense_matrix_t		= typename rom_data_t::dense_matrix_t;

  dense_matrix_t PhiTJPhi(romDim,romDim);
  PhiTJPhi.setZero();

  vector_t PhiTJPhiCol(romDim);
  PhiTJPhiCol.setZero();

  auto JPhiVecViewA = TestBasis.getMultiVectorView();
  //std::cout << " JPhiLocalLength in computeTestBasisTransposeTimesJacobianTimesBasis " <<  JPhiVecViewA.getLocalLength() << std::endl;


  computeTestBasisTransposeTimesJacobianTimesBasis_calc(appObj,y,mu,t,Phi,TestBasis,PhiTJPhi,PhiTJPhiCol,epsilon);
  return PhiTJPhi;
} 





template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeTestBasisTransposeTimesHessianTimesBasisTimesBasis_calc6(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const basis_t & TestBasis, hessian_t & BasisTransposeTimesHessianTimesBasisTimesBasis, vector_t & HCol,const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto ynorm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 

  //   - 45.0*im1jp3 - 405.0*im2jm1 + 81.0*im2jm2 - 9.0*im2jm3 + 405.0*im2jp1 - 81.0*im2jp2 + 9.0*im2jp3 + 45.0*im3jm1 - 9.0*im3jm2 + im3jm3 - 45.0*im3jp1 + 9.0*im3jp2 - 1.0*im3jp3 - 2025.0*ip1jm1 + 405.0*ip1jm2 - 45.0*ip1jm3 + 2025.0*ip1jp1 - 405.0*ip1jp2 + 45.0*ip1jp3 + 405.0*ip2jm1 - 81.0*ip2jm2 + 9.0*ip2jm3 - 405.0*ip2jp1 + 81.0*ip2jp2 - 9.0*ip2jp3 - 45.0*ip3jm1 + 9.0*ip3jm2 - 1.0*ip3jm3 + 45.0*ip3jp1 - 9.0*ip3jp2 + ip3jp3
  scalar_t scales[6];
  scales[0] =  -1./60.;
  scales[1] =   3./20.;
  scales[2] =  -3./4.;
  scales[3] =   3./4.;
  scales[4] =  -3./20.;
  scales[5] =   1./60.;

  int pert[6];
  pert[0] = -3;
  pert[1] = -2;
  pert[2] = -1;
  pert[3] = 1;
  pert[4] = 2;
  pert[5] = 3;

  for  (int i = 0; i < maxKForHessian; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      auto PhiCol_i = getCol<scalar_t>(Phi,i);
      auto PhiCol_j = getCol<scalar_t>(Phi,j);
      auto PhiColINorm = ::pressio::ops::norm2(PhiCol_i);
      auto PhiColJNorm = ::pressio::ops::norm2(PhiCol_j);
      scalar_t eps_i =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColINorm;
      scalar_t eps_j =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColJNorm;
      for (int k1=0; k1 < 6; k1++){
        for (int k2=0; k2 < 6; k2++){
          //first term
         ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i*pert[k1],PhiCol_j,eps_j*pert[k2]);
         ::pressio::ops::set_zero(ftmp);
          appObj.velocity(yp,t,ftmp);
          ::pressio::ops::update(workingVector,0.,ftmp,scales[k1]*scales[k2]);
        }
      }
      ::pressio::ops::product(::pressio::transpose(),1.,TestBasis,workingVector , 0., HCol);
      setSymmetricHessianCol(BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,i,j);
    }
  }
} 












template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeTestBasisTransposeTimesHessianTimesBasisTimesBasis_calc4(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const basis_t & TestBasis, hessian_t & BasisTransposeTimesHessianTimesBasisTimesBasis, vector_t & HCol,const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto ynorm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 
  for  (int i = 0; i < maxKForHessian; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      auto PhiCol_i = getCol<scalar_t>(Phi,i);
      auto PhiCol_j = getCol<scalar_t>(Phi,j);
      auto PhiColINorm = ::pressio::ops::norm2(PhiCol_i);
      auto PhiColJNorm = ::pressio::ops::norm2(PhiCol_j);
      scalar_t eps_i =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColINorm;
      scalar_t eps_j =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColJNorm;
      auto scale = 1./(144*eps_i*eps_j);

      //first term
      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i,PhiCol_j,-2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,0.,ftmp,8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,2.*eps_i,PhiCol_j,-eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-2*eps_i,PhiCol_j,eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-eps_i,PhiCol_j,2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,8.*scale);


      //second term
      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-eps_i,PhiCol_j,-2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-2.*eps_i,PhiCol_j,-1.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i,PhiCol_j,2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-8.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,2.*eps_i,PhiCol_j,eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-8.*scale);

      //third term
      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,2.*eps_i,PhiCol_j,-2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-2.*eps_i,PhiCol_j,2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-2.*eps_i,PhiCol_j,-2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,2.*eps_i,PhiCol_j,2.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);

 
      //fourth term
      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-1.*eps_i,PhiCol_j,-1.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,64.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i,PhiCol_j,eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,64.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,1.*eps_i,PhiCol_j,-1.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-64.*scale);

      ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,-1.*eps_i,PhiCol_j,1.*eps_j);
      ::pressio::ops::set_zero(ftmp);
      appObj.velocity(yp,t,ftmp);
      ::pressio::ops::update(workingVector,1.,ftmp,-64.*scale);

      ::pressio::ops::product(::pressio::transpose(),1.,TestBasis,workingVector , 0., HCol);
      setSymmetricHessianCol(BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,i,j);
    }
  }
} 

template <typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t, typename hessian_t, typename vector_t>
void computeTestBasisTransposeTimesHessianTimesBasisTimesBasis_calc(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi, const basis_t & TestBasis, hessian_t & BasisTransposeTimesHessianTimesBasisTimesBasis, vector_t & HCol,const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.updateScalarParameters(mu);
  appObj.velocity( y , t , f0);

  auto ynorm = ::pressio::ops::norm2(y);

  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 
  for  (int i = 0; i < maxKForHessian; i++){
    for (int j=0; j <= i; j++){
      ::pressio::ops::set_zero(workingVector);
      if (i == j){
        auto PhiCol = getCol<scalar_t>(Phi,i);
        auto PhiColNorm = ::pressio::ops::norm2(PhiCol);
        scalar_t eps =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColNorm;
        auto scale = 1./(12.*eps*eps);
        //yp = y + 2.*eps*PhiCol;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol,2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,-1.*scale);

        //yp = y + eps*PhiCol (or yp = yp - eps PhiCol);
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale);

    
        ::pressio::ops::update(workingVector,1.,f0,-30.*scale);
    
      
        //yp = y - eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-2.*eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,16.*scale);
       
        //yp = y - 2.*eps*PhiCol;
        ::pressio::ops::update(yp,1.,PhiCol,-eps);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
 
//        ::pressio::ops::scale(workingVector,1./(12.*(eps*eps)));
      }
      else{
        auto PhiCol_i = getCol<scalar_t>(Phi,i);
        auto PhiCol_j = getCol<scalar_t>(Phi,j);
        auto PhiColINorm = ::pressio::ops::norm2(PhiCol_i);
        auto PhiColJNorm = ::pressio::ops::norm2(PhiCol_j);
        scalar_t eps_i =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColINorm;
        scalar_t eps_j =  std::sqrt( (1. + ynorm)*epsilon ) / PhiColJNorm;
        auto scale = 0.25/(eps_i*eps_j);


        //yp = y + eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,0.,y,1.,PhiCol_i,eps_i,PhiCol_j,eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,0.,ftmp,1.*scale);

        //yp = y + eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);

        //yp = y - eps*PhiCol_i + eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_i,-2.*eps_i,PhiCol_j,2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,-1.*scale);
       
        //yp = y - eps*PhiCol_i - eps*PhiCol_j;
        ::pressio::ops::update(yp,1.,PhiCol_j,-2.*eps_j);
        ::pressio::ops::set_zero(ftmp);
        appObj.velocity(yp,t,ftmp);
        ::pressio::ops::update(workingVector,1.,ftmp,1.*scale);
 
      }
      ::pressio::ops::product(::pressio::transpose(),1., TestBasis,workingVector , 0., HCol);
      setSymmetricHessianCol(BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,i,j);
    }
  }
} 



  
template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesHessianTimesBasisTimesBasisKokkos(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,const basis_t & TestBasis, const int maxKForHessian, const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    =  typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesHessianTimesBasisTimesBasis("PhiTHPhiPhi", romDim,maxKForHessian,maxKForHessian);
  vector_t HCol("HCol",romDim);


  computeTestBasisTransposeTimesHessianTimesBasisTimesBasis_calc6(appObj,y,mu,t,Phi,TestBasis,BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,maxKForHessian,epsilon);
  return BasisTransposeTimesHessianTimesBasisTimesBasis; 
} 



template <typename rom_data_t, typename app_t, typename state_t, typename param_t , typename scalar_t, typename basis_t>
auto computeTestBasisTransposeTimesHessianTimesBasisTimesBasisEigen(const app_t & appObj,state_t & y,const param_t & mu,const scalar_t t,const basis_t & Phi,const basis_t & TestBasis, const int maxKForHessian,const scalar_t epsilon)
{
  auto romDim = ::pressio::ops::extent(Phi,1);
  auto fomDim = ::pressio::ops::extent(y,0);

  // Create temporary working vectors
  auto f0 = appObj.createVelocity();
  appObj.velocity( y , t , f0);
  auto ftmp = appObj.createVelocity();
  auto workingVector = appObj.createVelocity();

  auto yp = ::pressio::ops::clone(y); //pressio ops clone 

  using vector_t    = typename rom_data_t::vector_t;
  using dense_hessian_t    =  typename rom_data_t::dense_hessian_t;

  dense_hessian_t BasisTransposeTimesHessianTimesBasisTimesBasis(romDim);

  vector_t HCol(romDim);
  HCol.setZero();

  for (int i=0; i < romDim; i++){
    BasisTransposeTimesHessianTimesBasisTimesBasis[i].resize(maxKForHessian);
  }
  for (int i=0; i < romDim; i++){
    for (int j=0; j < maxKForHessian; j++){
      BasisTransposeTimesHessianTimesBasisTimesBasis[i][j].resize(maxKForHessian);
      for (int k=0; k < maxKForHessian; k++){
        BasisTransposeTimesHessianTimesBasisTimesBasis[i][j][k] = 0.;
      }
    }
  }

  computeTestBasisTransposeTimesHessianTimesBasisTimesBasis_calc6(appObj,y,mu,t,Phi,TestBasis,BasisTransposeTimesHessianTimesBasisTimesBasis,HCol,maxKForHessian,epsilon);
  return BasisTransposeTimesHessianTimesBasisTimesBasis; 
} 

}}}
#endif
