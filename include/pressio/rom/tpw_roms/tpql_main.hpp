#ifndef TPQL_MAIN_HPP_
#define TPQL_MAIN_HPP_

#include "tpql_operators.hpp"
#include "pressio/ops.hpp"

namespace pressio{ namespace rom{ namespace experimental{

template <typename app_t, typename basis_t, typename scalar_t>
class TpwqRom
{
public:
  template <typename list_of_states_t, typename list_of_params_t, 
            typename list_of_times_t, typename list_of_coords_t, typename fom_state_t> 
  TpwqRom(const app_t & appObj,const basis_t & Phi,const list_of_states_t & InputListOfStatesAtLinearizationPoints,
          const list_of_params_t & InputListOfParamsAtLinearizationPoints, const list_of_times_t & InputListOfTimesAtLinearizationPoints, 
          const list_of_coords_t & InputListOfCoordsAtLinearizationPoints, const fom_state_t & fomReferenceState, const scalar_t & velocityErrorTolerance,
          const scalar_t FiniteDifferenceStepSize) : Phi_(Phi) , appObj_(appObj){


  // Compute mean of coordinates
  int numberOfInputPoints = InputListOfCoordsAtLinearizationPoints.size(); 
  int coordinateDimension = InputListOfCoordsAtLinearizationPoints[0].size();

  std::vector<scalar_t> coordinateMeans(coordinateDimension);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension; j++){ 
      coordinateMeans[j] += InputListOfCoordsAtLinearizationPoints[i](j);
    }
  }

  for (int j=0; j < coordinateDimension; j++){ 
   coordinateMeans[j] /= float(numberOfInputPoints);
  }

  std::vector<scalar_t> distanceFromMean(numberOfInputPoints);
  for (int i = 0; i < numberOfInputPoints; i++){
    for (int j=0; j < coordinateDimension; j++){ 
      distanceFromMean[i] += std::pow( InputListOfCoordsAtLinearizationPoints[i](j) - coordinateMeans[j] ,2.);
    }
    distanceFromMean[i] = std::sqrt(distanceFromMean[i]);
  }

  std::vector<int> sortedIndices(numberOfInputPoints);
  for (int i=0;i<numberOfInputPoints;i++){sortedIndices[i] = i;};
  std::sort( sortedIndices.begin(),sortedIndices.end(), [&](int i,int j){return distanceFromMean[i]<distanceFromMean[j];} );
 

  auto workingStateVector = appObj.createVelocity(); 
  for (int i = 0; i < numberOfInputPoints; i++){
    // Construct state at current point
    int sortedIndex = sortedIndices[i];
    auto currentState = InputListOfStatesAtLinearizationPoints[sortedIndex];
    auto currentTime = InputListOfTimesAtLinearizationPoints[sortedIndex];
    auto currentParam = InputListOfParamsAtLinearizationPoints[sortedIndex];
    auto currentCoords = InputListOfCoordsAtLinearizationPoints[sortedIndex];

    ::pressio::ops::product(::pressio::nontranspose(), 1., Phi,currentState,0.,workingStateVector);
    ::pressio::ops::update(workingStateVector,1.,fomReferenceState,1.); 
    if (i == 0){
      // Add point to coordinate list
      ListOfStatesAtLinearizationPoints_.push_back( currentState );  
      ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
      ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
      ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  


      auto PhiTf = ::pressio::rom::experimental::computeBasisTransposeTimesVelocity(appObj, workingStateVector,currentParam,currentTime, Phi);
      ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);

      auto PhiTJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi);
      ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);

      auto PhiTHPhiPhi = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi);
      ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);

      auto PhiTJParams = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobian(appObj,workingStateVector,currentParam,currentTime, Phi);
      ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);

      auto PhiTHParams = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessian(appObj,workingStateVector,currentParam,currentTime, Phi);
      ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);

      auto PhiTHMixed = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessian(appObj,workingStateVector,currentParam,currentTime, Phi);
      ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);
      
    }
    else{

      auto relativeError = CheckGalerkinVelocity(currentState,workingStateVector,currentParam,currentTime,currentCoords);
        std::cout << " error = " << relativeError << ", passing point " << std::endl; 
      if (relativeError >= velocityErrorTolerance){
        std::cout << " error = " << relativeError << ", adding new point " << std::endl; 
        // Add point to coordinate list
        ListOfStatesAtLinearizationPoints_.push_back( currentState );  
        ListOfParamsAtLinearizationPoints_.push_back( currentParam );  
        ListOfTimesAtLinearizationPoints_.push_back( currentTime );  
        ListOfCoordsAtLinearizationPoints_.push_back( currentCoords );  
  

        auto PhiTf = ::pressio::rom::experimental::computeBasisTransposeTimesVelocity(appObj, workingStateVector,currentParam,currentTime, Phi);
        ListOfVelocitiesAtLinearizationPoints_.push_back(PhiTf);
  
        auto PhiTJPhi = ::pressio::rom::experimental::computeBasisTransposeTimesJacobianTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi);
        ListOfJacobiansAtLinearizationPoints_.push_back(PhiTJPhi);
  
        auto PhiTHPhiPhi = ::pressio::rom::experimental::computeBasisTransposeTimesHessianTimesBasisTimesBasis(appObj,workingStateVector,currentParam,currentTime, Phi);
        ListOfHessiansAtLinearizationPoints_.push_back(PhiTHPhiPhi);
  
        auto PhiTJParams = ::pressio::rom::experimental::computeBasisTransposeTimesParameterJacobian(appObj,workingStateVector,currentParam,currentTime, Phi);
        ListOfParamJacobiansAtLinearizationPoints_.push_back(PhiTJParams);
  
        auto PhiTHParams = ::pressio::rom::experimental::computeBasisTransposeTimesParameterHessian(appObj,workingStateVector,currentParam,currentTime, Phi);
        ListOfParamHessiansAtLinearizationPoints_.push_back(PhiTHParams);
  
        auto PhiTHMixed = ::pressio::rom::experimental::computeBasisTransposeTimesMixedHessian(appObj,workingStateVector,currentParam,currentTime, Phi);
        ListOfMixedHessiansAtLinearizationPoints_.push_back(PhiTHMixed);
  
      }
    }
  } 
 
  }
private:
  // types
  using rom_state_t = Eigen::Matrix<scalar_t,-1,1>;
  using params_t = Eigen::Matrix<scalar_t,-1,1>;
  using coords_t = Eigen::Matrix<scalar_t,-1,1>;

  using rom_velocity_t = Eigen::Matrix<scalar_t,-1,1>;
  using rom_jacobian_t = Eigen::Matrix<scalar_t,-1,-1>;
  using rom_hessian_t = std::vector< std::vector< std::vector< scalar_t> > >; 

  using list_of_states_t = std::vector<rom_state_t>;
  using list_of_params_t = std::vector<params_t>;
  using list_of_times_t = std::vector<scalar_t>;
  using list_of_coords_t = std::vector<coords_t>;

  const app_t & appObj_;
  const basis_t & Phi_;

  list_of_states_t ListOfStatesAtLinearizationPoints_{0};
  list_of_params_t ListOfParamsAtLinearizationPoints_{0}; 
  list_of_times_t  ListOfTimesAtLinearizationPoints_{0};
  list_of_coords_t  ListOfCoordsAtLinearizationPoints_{0};

  std::vector< rom_velocity_t > ListOfVelocitiesAtLinearizationPoints_{0};
  std::vector< rom_jacobian_t > ListOfJacobiansAtLinearizationPoints_{0};
  std::vector< rom_hessian_t  > ListOfHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_jacobian_t  > ListOfParamJacobiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfParamHessiansAtLinearizationPoints_{0}; 
  std::vector< rom_hessian_t  > ListOfMixedHessiansAtLinearizationPoints_{0}; 

  template <typename params_t, typename coords_t>
  rom_state_t GalerkinVelocity(const rom_state_t & xhat,const params_t & mu,const scalar_t t,const coords_t & PointCoords){
    int linPointIndx = FindClosestLinearizationPoint(PointCoords);
    auto VelocityAtLinearizationPoint = ListOfVelocitiesAtLinearizationPoints_[linPointIndx];
    auto JacobianAtLinearizationPoint = ListOfJacobiansAtLinearizationPoints_[linPointIndx];
    auto HessianAtLinearizationPoint =  ListOfHessiansAtLinearizationPoints_[linPointIndx];
    auto ParameterJacobianAtLinearizationPoint = ListOfParamJacobiansAtLinearizationPoints_[linPointIndx];
    auto ParameterHessianAtLinearizationPoint = ListOfParamHessiansAtLinearizationPoints_[linPointIndx];
    auto MixedHessianAtLinearizationPoint = ListOfMixedHessiansAtLinearizationPoints_[linPointIndx];
    auto ReducedStateAtLinearizationPoint = ListOfStatesAtLinearizationPoints_[linPointIndx];
    auto ParamAtLinearizationPoint = ListOfParamsAtLinearizationPoints_[linPointIndx];

    auto d_xHat = xhat - ReducedStateAtLinearizationPoint;
    auto d_mu = mu - ParamAtLinearizationPoint;
    auto f = VelocityAtLinearizationPoint + JacobianAtLinearizationPoint * d_xHat + 
             ParameterJacobianAtLinearizationPoint * d_mu;
    tensorMultiply(HessianAtLinearizationPoint, d_xHat,d_xHat,0.5,f); 
    tensorMultiply(ParameterHessianAtLinearizationPoint, d_mu,d_mu,0.5,f); 
    tensorMultiply(MixedHessianAtLinearizationPoint, d_xHat,d_mu,1.,f); 
    return f;
  }

 
  void tensorMultiply(rom_hessian_t H, rom_state_t a1, rom_state_t a2,scalar_t alpha, rom_velocity_t f){
    auto d1 = H.size();
    auto d2 = H[0].size();
    auto d3 = H[0][0].size();
    for (int i = 0; i < d1; i++){
      for (int j = 0; j < d2; j++){
        for (int k = 0; k < d3; k++){
          f[i] += alpha*H[i][j][k] * a1(j) * a2(k);
        }
      }
    } 
  }

  template <typename fom_state_t, typename params_t, typename coords_t>
  double CheckGalerkinVelocity(const rom_state_t & xhat,const fom_state_t & x, const params_t & mu,const scalar_t t,const coords_t & PointCoords)
  {
    auto fTpwq = GalerkinVelocity(xhat,mu,t,PointCoords);
    auto fTrue = appObj_.createVelocity();
    //appObj.updateParams(mu);
    appObj_.velocity(x,t,fTrue);

    rom_velocity_t fGalerkin( xhat.size());
    ::pressio::ops::product(::pressio::transpose(),1., Phi_,fTrue, 0.,fGalerkin);    
    auto error = fTpwq - fGalerkin;
    auto errorNorm = error.norm() / fGalerkin.norm(); 
    return errorNorm; 
  }

  template<typename coords_t>
  int FindClosestLinearizationPoint(coords_t PointCoords){

    int numberOfAddedLinearizationPoints = ListOfCoordsAtLinearizationPoints_.size();
    int coordinateDimension = 2;
    std::vector<scalar_t> distance(numberOfAddedLinearizationPoints);
    std::cout << "Current coords are " << PointCoords << std::endl;

    for (int i = 0; i < numberOfAddedLinearizationPoints; i++){
      for (int j=0; j < coordinateDimension; j++){
        distance[i] += std::pow( ListOfCoordsAtLinearizationPoints_[i](j) - PointCoords[j] ,2.);
      }
      distance[i] = std::sqrt(distance[i]);
    }
    int linPointIndex = std::min_element(distance.begin(), distance.end()) - distance.begin();
    std::cout << "Closest linearization point is " << linPointIndex << std::endl;
    return linPointIndex;
  }

  //template <typename params_t, typename coords_t>
  //void GalerkinVelocity(const rom_state_t & xhat,const params_t & mu,const scalar_t t,const coords_t & PointCoords){
  //}
};
}}}
#endif
