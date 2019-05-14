
#ifndef ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_
#define ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_

#include "../../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace rompp{ namespace rom{ namespace impl{

//#ifdef HAVE_TRILINOS
// template <
//   ode::ImplicitEnum odeMethod,
//   typename lspg_matrix_type,
//   typename scalar_type,
//   typename decoder_jac_type,
//   ::rompp::mpl::enable_if_t<
//     core::meta::is_multi_vector_wrapper_tpetra_block<lspg_matrix_type>::value and
//     core::meta::is_multi_vector_wrapper_tpetra_block<decoder_jac_type>::value
//     > * = nullptr
//   >
// void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi stands for J * phi
// 			    scalar_type	dt,
// 			    const decoder_jac_type & phi){

//   constexpr auto one = ::rompp::core::constants::one<scalar_type>();
//   constexpr auto negOne = ::rompp::core::constants::negOne<scalar_type>();
//   auto coeff = negOne*dt;
//   if (odeMethod == ::rompp::ode::ImplicitEnum::BDF2)
//     coeff *= ::rompp::ode::coeffs::bdf2<scalar_type>::c3_;

//   jphi.data()->update( one, *phi.data(), coeff);
// }
// #endif


template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_eigen<lspg_matrix_type>::value and
    core::meta::is_multi_vector_wrapper_eigen<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  //loop over elements of jphi
  for (auto i=0; i<jphi.length(); i++){
    for (auto j=0; j<jphi.numVectors(); j++)
      jphi(i,j) = phi(i,j) - prefactor*dt*jphi(i,j);
  }
}



/* when we have hyper-reduction and we need to calculate the
 * time-discrete residual, yn (i.e. the current state) and
 * ynm (i.e. the states at prev steps) might have
 * different entries/maps than R (i.e. the spatial residual)
 * so we need to update R only picking elements corresponding
 * elements in yn and ynm.
 * In other words, in the functions below we do NOT assume that
 * yn and ynm have the same distribution/sizes of R but we DO assume
 * that yn and ynm contain at least all the elements of R (and more).
 * THis is because to evaluate spatial residual one typically
 * needs neighbors values, so in the hyper-red scenario the
 * states contain more DOFs than the spatial residual.
 * for instance, imagine the i-th cell in a uniform grid with a
 * centered 2nd-order FD operator.
 * To compute the spatial residual at the i-th cell, one needs the
 * cell value y_i, as well as the left, y_i-1, and right, y_i+1,
 * values of the state. So for a single DOF of the residual, one
 * needs 3 DOFs of the state. Extrapolating to all cells, then we
 * clearly see that we need more dofs for the states than for the residual.
 *
 * Below, we account for this potential scenario and the code below
 * obviously also works for the case where the vectors have same
 * the same distributions.
*/

#ifdef HAVE_TRILINOS

template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_epetra<lspg_matrix_type>::value and
    core::meta::is_multi_vector_wrapper_epetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi stands for J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  // integral type of the global indices
  using GO_t = typename core::details::traits<lspg_matrix_type>::global_ordinal_t;

  //get row map of phi
  const auto & phi_map = phi.getDataMap();
  // get my global elements
  std::vector<GO_t> gIDphi( phi.localLength() );
  phi_map.MyGlobalElements( gIDphi.data() );

  // get map of jphi
  const auto & jphi_map = jphi.getDataMap();
  // get global elements
  std::vector<GO_t> gIDjphi( jphi.localLength() );
  jphi_map.MyGlobalElements( gIDjphi.data() );

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  //loop over elements of jphi
  for (auto i=0; i<jphi.localLength(); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    auto lid = phi_map.LID(gIDjphi[i]);
    for (auto j=0; j<jphi.globalNumVectors(); j++)
      jphi(i,j) = phi(lid,j) - prefactor*dt*jphi(i,j);
  }
}

template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_tpetra<lspg_matrix_type>::value and
    core::meta::is_multi_vector_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  //get row map of phi
  const auto phi_map = phi.getMap();
  // get my global elements
  const auto gIDphi = phi_map->getMyGlobalIndices();

  // get map of jphi
  const auto jphi_map = jphi.getMap();
  // get global elements
  const auto gIDjphi = jphi_map->getMyGlobalIndices();

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  auto jphi2dView = jphi.get2dViewNonConst();
  auto phi2dView = phi.get2dView();

  //loop over elements of jphi
  for (size_t i=0; i<jphi.getLocalLength(); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    const auto lid = phi_map->getLocalElement(gIDjphi[i]);
    for (size_t j=0; j<jphi.getNumVectors(); j++)
      jphi2dView[j][i] = phi2dView[j][lid] - prefactor*dt*jphi2dView[j][i];
  }
}

template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra<lspg_matrix_type>::value and
    core::meta::is_multi_vector_wrapper_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi,
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  //time_discrete_jacobian<odeMethod>(*jphi.data(), dt, *phi.data());
}


template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra_block<lspg_matrix_type>::value and
    core::meta::is_multi_vector_wrapper_tpetra_block<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi,
			    scalar_type	dt,
			    const decoder_jac_type & phi){
  // auto jphi_mvv = jphi.data()->getMultiVectorView();
  // auto phi_mvv  = phi.data()->getMultiVectorView();
  // time_discrete_jacobian<odeMethod>(jphi_mvv, dt, phi_mvv);
}

#endif

}}}//end namespace rompp::rom::impl
#endif
