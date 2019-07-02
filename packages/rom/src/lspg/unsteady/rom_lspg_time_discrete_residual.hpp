
#ifndef ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_
#define ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_

#include "../../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace pressio{ namespace rom{ namespace impl{


template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  typename ud_ops,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::Euler 
#ifdef HAVE_PYBIND11 
    and mpl::not_same< ud_ops, pybind11::object>::value
#endif    
   > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt,
			    const ud_ops * udOps){

  udOps->time_discrete_euler(*R.data(), *yn.data(),
			     *ynm[0].data(), dt);
}

#ifdef HAVE_PYBIND11
template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  typename ud_ops,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::Euler and
    mpl::is_same< ud_ops, pybind11::object>::value
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt,
			    const ud_ops & udOps){

  udOps.attr("time_discrete_euler")(R, yn, ynm[0], dt);
}
#endif    




/*
 * for EIGEN
*/
template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_eigen<state_type>::value == true and
    method == ::pressio::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)
  *R.data() = *yn.data() - *ynm[0].data() - dt * (*R.data());

}

template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_eigen<state_type>::value == true and
    method == ::pressio::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  using namespace ::pressio::ode::coeffs;

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
template <::pressio::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_epetra;

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    R = yn[lid] - ynm[0][lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::pressio::ode::coeffs;

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
  ::pressio::ode::ImplicitEnum odeMethod,
  int numStates,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_epetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,numStates> & ynm,
			    state_type & R,
			    scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)

  // the integral type of the global indices
  using GO_t = typename containers::details::traits<state_type>::global_ordinal_t;

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




/*************************************
	             tpetra
*************************************/
template <::pressio::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_tpetra;

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    //R[i] = yn[lid]- ynm[0][lid] - dt*R[i];
    R = (yn.getData())[lid]
      - (ynm[0].getData())[lid] - dt * R;
  }

  template <typename T, typename state_type>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const state_type & ynm0){
    R = (yn.getData())[lid]
      - (ynm0.getData())[lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type & yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::pressio::ode::coeffs;

    // R[i] = yn[lid]
    //   - bdf2<scalar_type>::c1_*ynm[1][lid]
    //   - bdf2<scalar_type>::c2_*ynm[0][lid]
    //   - bdf2<scalar_type>::c3_*dt*R[i];
    R = (yn.getData())[lid]
      - bdf2<T>::c1_*(ynm[1].getData())[lid]
      + bdf2<T>::c2_*(ynm[0].getData())[lid]
      - bdf2<T>::c3_*dt*R;
  }

  template <typename T, typename state_type>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type& yn,
			 const state_type & ynm0,
			 const state_type & ynm1){
    using namespace ::pressio::ode::coeffs;
    R = (yn.getData())[lid]
      - bdf2<T>::c1_*(ynm1.getData())[lid]
      + bdf2<T>::c2_*(ynm0.getData())[lid]
      - bdf2<T>::c3_*dt*R;
  }
};


template<
  ::pressio::ode::ImplicitEnum odeMethod,
  typename state_type,
  typename scalar_type,
  typename ... Args,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_tpetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual_tpetra_impl(const state_type & yn,
			    state_type & R,
			    scalar_type dt,
			    Args && ... args)/*
			    const state_type & ynm0,
			    const state_type & ynm1,)*/
{
  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)

  // get map of yn (ynm has for sure the same map as yn)
  const auto y_map = yn.getMap();
  // get my global elements
  const auto gIDy = y_map->getMyGlobalIndices();

  // get map of R
  const auto R_map = R.getMap();
  // get global elements
  const auto gIDr = R_map->getMyGlobalIndices();

  //loop over elements of R
  for (size_t i=0; i<R.getLocalLength(); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    const auto lid = y_map->getLocalElement(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_tpetra<
      odeMethod
      >::template evaluate<
      scalar_type, state_type
      >(dt, R.getDataNonConst()[i], lid, yn,
	std::forward<Args>(args)...);//ynm0, ynm1);
  }
}


template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    odeMethod == ::pressio::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  time_discrete_residual_tpetra_impl<odeMethod>
    (*yn.data(),
     *R.data(), dt,
     *ynm[0].data());

}

template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    odeMethod == ::pressio::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  time_discrete_residual_tpetra_impl<odeMethod>
    (*yn.data(),
     *R.data(), dt,
     *ynm[0].data(),
     *ynm[1].data());
}




/*************************************
            tpetra block
*************************************/
template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    odeMethod == ::pressio::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  auto yn_vv = const_cast<state_type &>(yn).data()->getVectorView();
  auto ynm0_vv = const_cast<state_type &>(ynm[0]).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<odeMethod>(
      yn_vv, R_vv, dt, ynm0_vv);
}

template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    odeMethod == ::pressio::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  auto yn_vv = const_cast<state_type &>(yn).data()->getVectorView();
  auto ynm0_vv = const_cast<state_type &>(ynm[0]).data()->getVectorView();
  auto ynm1_vv = const_cast<state_type &>(ynm[1]).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<odeMethod>(
    yn_vv, R_vv, dt, ynm0_vv, ynm1_vv);
}

#endif

}}}//end namespace pressio::rom::impl
#endif
