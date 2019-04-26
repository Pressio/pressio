
#ifndef ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_
#define ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_

#include "../../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace rompp{ namespace rom{ namespace impl{

#ifdef HAVE_TRILINOS
/* For block:Tpetra, for time being, just implement Euler also assuming that
 * there are no gaps in the residual and states, so basically
 * assuming that states and residual have same distributions.
 * Tis has to be fixed later similarly to what we do for regular-tpetra and epetra*/
template<
  ::rompp::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    odeMethod == ::rompp::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(state_type & yn, // there is a const missing
			    std::array<state_type,n> & ynm,// there is a const missing
			    state_type & R,
			    scalar_type dt){
  constexpr auto one = ::rompp::core::constants::one<scalar_type>();
  constexpr auto negOne = ::rompp::core::constants::negOne<scalar_type>();
  const scalar_type negDt = -dt;
  // R.data()->update( 1.0, *yn.data(), negDt);
  // R.data()->update( -1.0, *ynm[0].data(), one);
  ::rompp::core::ops::do_update(R, negDt, yn, one, ynm[0], negOne);
}
#endif

template<
  ::rompp::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_eigen<state_type>::value == true and
    method == ::rompp::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)
  *R.data() = *yn.data() - *ynm[0].data() - dt * (*R.data());
}


template<
  ::rompp::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_eigen<state_type>::value == true and
    method == ::rompp::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  using namespace ::rompp::ode::coeffs;

  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)
  R = yn
    - bdf2<scalar_type>::c1_*ynm[1]
    + bdf2<scalar_type>::c2_*ynm[0]
    - bdf2<scalar_type>::c3_*dt*R;
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
template <::rompp::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_epetra;

template <>
struct time_discrete_single_entry_epetra<::rompp::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    //R[i] = yn[lid]- ynm[0][lid] - dt*R[i];
    R = yn[lid] - ynm[0][lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_epetra<::rompp::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::rompp::ode::coeffs;

    // R[i] = yn[lid]
    //   - bdf2<scalar_type>::c1_*ynm[1][lid]
    //   - bdf2<scalar_type>::c2_*ynm[0][lid]
    //   - bdf2<scalar_type>::c3_*dt*R[i];
    R = yn[lid]
      - bdf2<T>::c1_*ynm[1][lid]
      + bdf2<T>::c2_*ynm[0][lid]
      - bdf2<T>::c3_*dt*R;
  }
};


template<
  ::rompp::ode::ImplicitEnum odeMethod,
  int numStates,
  typename state_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_epetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,numStates> & ynm,
			    state_type & R,
			    scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)

  // the integral type of the global indices
  using GO_t = typename core::details::traits<state_type>::global_ordinal_t;

  // get map of yn (ynm has for sure the same map as yn)
  const auto & y_map = yn.getDataMap();
  // get my global elements
  std::vector<GO_t> gIDy( yn.localSize() );
  y_map.MyGlobalElements( gIDy.data() );

  // get map of R
  const auto & R_map = R.getDataMap();
  // get global elements
  std::vector<GO_t> gIDr( R.localSize() );
  R_map.MyGlobalElements( gIDr.data() );

  //loop over elements of R
  for (auto i=0; i<R.localSize(); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    auto lid = y_map.LID(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_epetra<
      odeMethod
      >::template evaluate<
      scalar_type, state_type, numStates>(dt, R[i], lid, yn, ynm);
  }
}



template <::rompp::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_tpetra;

template <>
struct time_discrete_single_entry_tpetra<::rompp::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    //R[i] = yn[lid]- ynm[0][lid] - dt*R[i];
    R = (yn.data()->getData())[lid]
      - (ynm[0].data()->getData())[lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_tpetra<::rompp::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::rompp::ode::coeffs;

    // R[i] = yn[lid]
    //   - bdf2<scalar_type>::c1_*ynm[1][lid]
    //   - bdf2<scalar_type>::c2_*ynm[0][lid]
    //   - bdf2<scalar_type>::c3_*dt*R[i];
    R = (yn.data()->getData())[lid]
      - bdf2<T>::c1_*(ynm[1].data()->getData())[lid]
      + bdf2<T>::c2_*(ynm[0].data()->getData())[lid]
      - bdf2<T>::c3_*dt*R;
  }
};

template<
  ::rompp::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_tpetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)

  // get map of yn (ynm has for sure the same map as yn)
  const auto & y_map = yn.getDataMap();
  // get my global elements
  auto gIDy = y_map.getMyGlobalIndices();

  // get map of R
  const auto & R_map = R.getDataMap();
  // get global elements
  auto gIDr = R_map.getMyGlobalIndices();

  //loop over elements of R
  for (auto i=0; i<R.localSize(); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    auto lid = y_map.getLocalElement(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_tpetra<
      odeMethod
      >::template evaluate<
      scalar_type, state_type,
      n>(dt, (R.data()->getDataNonConst())[i], lid, yn, ynm);
  }
}

#endif

}}}//end namespace rompp::rom::impl
#endif
