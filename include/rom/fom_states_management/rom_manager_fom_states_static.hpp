/*
//@HEADER
// ************************************************************************
//
// rom_manager_fom_states_static.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ROM_FOM_STATES_MANAGEMENT_ROM_MANAGER_FOM_STATES_STATIC_HPP_
#define ROM_FOM_STATES_MANAGEMENT_ROM_MANAGER_FOM_STATES_STATIC_HPP_

namespace pressio{ namespace rom{

template <
  typename tag,
  typename fom_state_type,
  typename reconstuctor_type,
  typename ud_ops_t,
  std::size_t numstates = 1
  >
class ManagerFomStates;

// *****************************************
// partial specialize for steady
// *****************************************
template <
  typename fom_state_type,
  typename reconstuctor_type,
  typename ud_ops_t
  >
class ManagerFomStates<
  ::pressio::rom::Steady,
  fom_state_type, reconstuctor_type, ud_ops_t, 1
  >
{
public:
  using data_type =
    ::pressio::containers::IndexableStaticCollection<fom_state_type, 1>;
  using value_type = fom_state_type;

  ManagerFomStates() = delete;
  ManagerFomStates(const ManagerFomStates &) = default;
  ManagerFomStates & operator=(const ManagerFomStates &) = delete;
  ManagerFomStates(ManagerFomStates &&) = default;
  ManagerFomStates & operator=(ManagerFomStates &&) = delete;
  ~ManagerFomStates() = default;

  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

  template<
    typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<!std::is_void<_ud_ops_t>::value, int> = 0
    >
  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState,
		   const _ud_ops_t & udOps)
    : udOps_(udOps),
      fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

public:
  static constexpr std::size_t size(){ return data_type::size(); }

  fom_state_type const & currentFomState() const {return data_(0);}

  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romStateIn)
  {
    fomStateReconstrObj_(romStateIn, data_(0));
  }

private:
  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    ::pressio::ops::set_zero(data_(0));
  }

  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    udOps_.get().set_zero(*data_(0).data());
  }

private:
  typename std::conditional<
  std::is_void<ud_ops_t>::value,
  ::pressio::utils::impl::empty,
  std::reference_wrapper<const ud_ops_t>
  >::type udOps_;

  std::reference_wrapper<const reconstuctor_type> fomStateReconstrObj_;
  data_type data_;
};


// *****************************************
// partial specialize for explicit stepping
// *****************************************
template <
  typename fom_state_type,
  typename reconstuctor_type,
  typename ud_ops_t,
  std::size_t numstates
  >
class ManagerFomStates<
  ::pressio::rom::UnsteadyExplicit, fom_state_type, reconstuctor_type,
  ud_ops_t, numstates
  >
{
  static_assert
  (numstates>=1, "ManagerFomStates cannot be empty.");

public:
  using data_type  = ::pressio::containers::IndexableStaticCollection<fom_state_type,
								      numstates>;
  using value_type = fom_state_type;

  ManagerFomStates() = delete;
  ManagerFomStates(const ManagerFomStates &) = default;
  ManagerFomStates & operator=(const ManagerFomStates &) = delete;
  ManagerFomStates(ManagerFomStates &&) = default;
  ManagerFomStates & operator=(ManagerFomStates &&) = delete;
  ~ManagerFomStates() = default;

  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

  template<
    typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<!std::is_void<_ud_ops_t>::value, int> = 0
    >
  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState,
		   const _ud_ops_t & udOps)
    : udOps_(udOps),
      fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

public:
  static constexpr std::size_t size(){ return data_type::size(); }

  // for explicit time stepping it makese sense to have indexing
  // using "n, n-1, n-2", etc

  //---------------------------------
  // ** methods to const view data **
  //---------------------------------
  // n
  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=1, fom_state_type const &>
  fomStateAt(::pressio::ode::n) const {return data_(0);}

  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=1, fom_state_type const &>
  operator()(::pressio::ode::n) const {return data_(0);}

  // n-1
  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=2, fom_state_type const &>
  fomStateAt(::pressio::ode::nMinusOne) const {return data_(1);}

  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=2, fom_state_type const &>
  operator()(::pressio::ode::nMinusOne) const {return data_(1);}

  //-----------------------------
  // ** reconstruction methods **
  //-----------------------------
  /* for explicit methods, the current FOM state is the n-th state */
  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romStateIn)
  {
    this->reconstructAt(romStateIn, ::pressio::ode::n());
  }

  // n
  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=1 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::n)
  {
    fomStateReconstrObj_(romStateIn, data_(0));
  }

  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=1 >
  reconstructFomStateAt(const rom_state_t & romStateIn,
			::pressio::ode::n tag)
  {
    reconstructAt(romStateIn, tag);
  }

  // n-1
  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=2 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::n)
  {
    fomStateReconstrObj_(romStateIn, data_(1));
  }

  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=2 >
  reconstructFomStateAt(const rom_state_t & romStateIn,
			::pressio::ode::n tag)
  {
    reconstructAt(romStateIn, tag);
  }

private:
  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      ::pressio::ops::set_zero(data_(i));
  }

  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      udOps_.get().set_zero(*data_(i).data());
  }

private:
  typename std::conditional<
  std::is_void<ud_ops_t>::value,
  ::pressio::utils::impl::empty,
  std::reference_wrapper<const ud_ops_t>
  >::type udOps_;

  std::reference_wrapper<const reconstuctor_type> fomStateReconstrObj_;
  data_type data_;
};


// *****************************************
// partial specialize for implicit stepping
// *****************************************
template <
  typename fom_state_type,
  typename reconstuctor_type,
  typename ud_ops_t,
  std::size_t numstates
  >
class ManagerFomStates<
  ::pressio::rom::UnsteadyImplicit, fom_state_type, reconstuctor_type,
  ud_ops_t, numstates
  >
{
  static_assert
  (numstates>=1, "ManagerFomStates cannot be empty.");

public:
  using data_type  = ::pressio::containers::IndexableStaticCollection<fom_state_type,
								      numstates>;
  using value_type = fom_state_type;

  ManagerFomStates() = delete;
  ManagerFomStates(const ManagerFomStates &) = default;
  ManagerFomStates & operator=(const ManagerFomStates &) = delete;
  ManagerFomStates(ManagerFomStates &&) = default;
  ManagerFomStates & operator=(ManagerFomStates &&) = delete;
  ~ManagerFomStates() = default;

  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

  template<
    typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<!std::is_void<_ud_ops_t>::value, int> = 0
    >
  ManagerFomStates(const reconstuctor_type & fomStateReconstr,
		   const fom_state_type & fomState,
		   const _ud_ops_t & udOps)
    : udOps_(udOps),
      fomStateReconstrObj_(fomStateReconstr),
      data_(fomState)
  {
    this->setContainersToZero();
  }

public:
  static constexpr std::size_t size(){ return data_type::size(); }

  // for implicit time stepping it makese sense to
  // index using "n+1, n, n-1, n-2", etc

  // ** methods to extract const ref to data **
  // n+1
  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=1, fom_state_type const &>
  fomStateAt(::pressio::ode::nPlusOne) const {return data_(0);}

  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=1, fom_state_type const &>
  operator()(::pressio::ode::nPlusOne) const {return data_(0);}

  // n
  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=2, fom_state_type const &>
  fomStateAt(::pressio::ode::n) const {return data_(1);}

  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=2, fom_state_type const &>
  operator()(::pressio::ode::n) const {return data_(1);}

  // n-1
  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=3, fom_state_type const &>
  fomStateAt(::pressio::ode::nMinusOne) const {return data_(2);}

  template <std::size_t _numstates = numstates>
  mpl::enable_if_t<_numstates>=3, fom_state_type const &>
  operator()(::pressio::ode::nMinusOne) const {return data_(2);}

  //----------------------------------------
  // ** methods to reconstruct fom state **
  //----------------------------------------
  /* for implicit stepping, the current FOM state is the n+1-th state */
  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romStateIn)
  {
    this->reconstructAt(romStateIn, ::pressio::ode::nPlusOne());
  }

  // n+1
  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=1 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::nPlusOne)
  {
    fomStateReconstrObj_.get()(romStateIn, data_(0));
  }

  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=1 >
  reconstructFomStateAt(const rom_state_t & romStateIn,
			::pressio::ode::nPlusOne tag)
  {
    reconstructAt(romStateIn, tag);
  }

  // n
  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=2 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::n)
  {
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates>=2 >
  reconstructFomStateAt(const rom_state_t & romStateIn,
			::pressio::ode::n tag)
  {
    reconstructAt(romStateIn, tag);
  }

  //-----------------------------
  // reconstruct with update:
  //-----------------------------
  // reconstructs at point and shifts back existing FOM states
  // so that stencil is updating properly
  // we do this from n since n+1 is handled differenetly

  // n==2 we have y_n+1, y_n
  template <typename rom_state_t, std::size_t _numstates = numstates>
  mpl::enable_if_t< _numstates==2 >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /* when n == 2, it means I only have n+1 and n
     * so to reconstruct at n, I can simply
     * overwrite the data in data_(1) */
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  /* when n == 3, we have y_n+1, y_n, y_n-1 */
  template <
    typename rom_state_t,
    typename _ud_ops_t = ud_ops_t,
    std::size_t _numstates = numstates
    >
  mpl::enable_if_t< _numstates==3 and std::is_void<_ud_ops_t>::value >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /*
     * copy y_n into y_n-1
     * then reconstruct y_n
     */
    ::pressio::ops::deep_copy(data_(2), data_(1));
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  template <
    typename rom_state_t,
    typename _ud_ops_t = ud_ops_t,
    std::size_t _numstates = numstates
    >
  mpl::enable_if_t< _numstates==3 and !std::is_void<_ud_ops_t>::value >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    udOps_.get().deep_copy(*data_(2).data(), *data_(1).data());
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  /* when n == 4, we have y_n+1, y_n, y_n-1, y_n-2 */
  template <
    typename rom_state_t,
    typename _ud_ops_t = ud_ops_t,
    std::size_t _numstates = numstates
    >
  mpl::enable_if_t< _numstates==4 and std::is_void<_ud_ops_t>::value >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /*
     * copy y_n-1 into y_n-2
     * copy y_n into y_n-1
     * then reconstruct y_n */
    ::pressio::ops::deep_copy(data_(3), data_(2));
    ::pressio::ops::deep_copy(data_(2), data_(1));
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  template <
    typename rom_state_t,
    typename _ud_ops_t = ud_ops_t,
    std::size_t _numstates = numstates
    >
  mpl::enable_if_t< _numstates==4 and !std::is_void<_ud_ops_t>::value >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    udOps_.get().deep_copy(*data_(3).data(), *data_(2).data());
    udOps_.get().deep_copy(*data_(2).data(), *data_(1).data());
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

private:
  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      ::pressio::ops::set_zero(data_(i));
  }

  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void< _ud_ops_t>::value >
  setContainersToZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      udOps_.get().set_zero(*data_(i).data());
  }

private:
  typename std::conditional<
  std::is_void<ud_ops_t>::value,
  ::pressio::utils::impl::empty,
  std::reference_wrapper<const ud_ops_t>
  >::type udOps_;

  std::reference_wrapper<const reconstuctor_type> fomStateReconstrObj_;
  data_type data_;
};

}}//end namespace pressio::rom















// template <
//   std::size_t n,
//   typename fom_state_type,
//   typename reconstuctor_type,
//   typename ud_ops_t
//   >
// class ManagerFomStatesStatic
// {
//   static_assert
//   ( ::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
//     "ManagerFomStatesStatic: the fom_state_type must be a pressio wrapper.");

//   static_assert
//   (n>=1, "ManagerFomStatesStatic cannot be empty.");

// public:
//   using data_type  = ::pressio::containers::IndexableStaticCollection<fom_state_type, n>;
//   using value_type = fom_state_type;

//   ManagerFomStatesStatic() = delete;
//   ManagerFomStatesStatic(const ManagerFomStatesStatic &) = default;
//   ManagerFomStatesStatic & operator=(const ManagerFomStatesStatic &) = delete;
//   ManagerFomStatesStatic(ManagerFomStatesStatic &&) = default;
//   ManagerFomStatesStatic & operator=(ManagerFomStatesStatic &&) = delete;
//   ~ManagerFomStatesStatic() = default;

//   ManagerFomStatesStatic(const reconstuctor_type & fomStateReconstr,
// 			 const fom_state_type & fomState)
//     : fomStateReconstrObj_(fomStateReconstr),
//       data_(fomState)
//   {
//     PRESSIOLOG_DEBUG
//       ("cnstr: allocating n = {} fomStates with extent = {}",
//        n, fomState.extent(0));

//     this->resetContainersToZero();
//   }

//   template<
//     typename _ud_ops_t = ud_ops_t,
//     mpl::enable_if_t<!std::is_void<_ud_ops_t>::value, int> = 0
//     >
//   ManagerFomStatesStatic(const reconstuctor_type & fomStateReconstr,
// 			 const fom_state_type & fomState,
// 			 const _ud_ops_t & udOps)
//     : udOps_(udOps),
//       fomStateReconstrObj_(fomStateReconstr),
//       data_(fomState)
//   {
//     PRESSIOLOG_DEBUG
//       ("cnstr: allocating n = {} fomStates with extent = {}",
//        n, fomState.extent(0));

//     this->resetContainersToZero();
//   }

// public:
//   static constexpr std::size_t size() {
//     return data_type::size();
//   }

//   const fom_state_type & currentFomStateCRef() const{
//     static_assert( n>=1, "Cannot call currentFomStateCRef if n < 1");
//     return data_(0);
//   }

//   const fom_state_type & fomStatePrevStepCRef() const{
//     static_assert( n>=2, "Cannot call fomStatePrevStepCRef if n < 2");
//     return data_(1);
//   }

//   const fom_state_type & fomStatePrevPrevStepCRef() const{
//     static_assert( n>=3, "Cannot call fomStatePrevPrevStepCRef if n < 3");
//     return data_(2);
//   }

//   /* this method reconstructs the current FOM state */
//   template <typename rom_state_t>
//   void reconstructCurrentFomState(const rom_state_t & romStateIn)
//   {
// #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
//     auto timer = Teuchos::TimeMonitor::getStackedTimer();
//     timer->start("reconstruct fom state");
// #endif

//     static_assert( n>=1, "Cannot call reconstructCurrentFomState if n < 1");
//     fomStateReconstrObj_.get()(romStateIn, data_(0));

// #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
//     timer->stop("reconstruct fom state");
// #endif
//   }

//   // overload the left shift operator to use when we need to
//   // reconstruct the FOM state at n-1 and shift all previous ones

//   /* when n == 1, disable the << operator since that is meant
//    * to use for time-dep problems when we have to store the states
//    * for the stepper stencil */

//   /* when n == 2, it means I only have current state and previous one
//    * so when I need to reconstruct the previous state, I can simply
//    * overwrite the data in data_(1) */
//   template < typename rom_state_t, std::size_t _n = n>
//   ::pressio::mpl::enable_if_t<_n == 2>
//   operator << (const rom_state_t & romStateIn)
//   {
//     // reconstrct the FOM state at n-1
//     fomStateReconstrObj_.get()(romStateIn, data_(1));
//   }

//   /* when n >= 3, we need to deep copy data to
//    * such that y_t-2 goes into y_t-3,
//    * and y_t-1 goes into y_t-2, etc.
//    * and then finally we overwrite data_(0) */
//   template <typename _ud_ops_t = ud_ops_t, typename rom_state_t, std::size_t _n = n >
//   ::pressio::mpl::enable_if_t< _n >= 3 and std::is_void<_ud_ops_t>::value >
//   operator << (const rom_state_t & romStateIn)
//   {
//     // copy all states back, such that y_t-2 goes into y_t-3,
//     // and y_t-1 goes into y_t-2, etc. so that y_t-1 is free to overwrite
//     for (std::size_t i=n-2; i>=1; --i){
//       const auto & src  = data_(i);
//       auto & dest = data_(i+1);
//       ::pressio::ops::deep_copy(dest, src);
//     }
//     // then, reconstrct the FOM state at t-1
//     fomStateReconstrObj_.get()(romStateIn, data_(1));
//   }

//   template <typename _ud_ops_t = ud_ops_t, typename rom_state_t, std::size_t _n = n >
//     ::pressio::mpl::enable_if_t< _n >= 3 and !std::is_void<_ud_ops_t>::value >
//   operator << (const rom_state_t & romStateIn)
//   {
//     for (std::size_t i=n-2; i>=1; --i){
//       const auto & src  = data_(i);
//       auto & dest = data_(i+1);
//       udOps_.get().deep_copy(*dest.data(), *src.data());
//     }
//     // then, reconstrct the FOM state at t-1
//     fomStateReconstrObj_.get()(romStateIn, data_(1));
//   }

// private:
//   template <typename _ud_ops_t = ud_ops_t>
//   mpl::enable_if_t< std::is_void< _ud_ops_t>::value >
//   resetContainersToZero(){
//     for (std::size_t i=0; i<n; i++)
//       ::pressio::ops::set_zero(data_(i));
//   }

//   template <typename _ud_ops_t = ud_ops_t>
//   mpl::enable_if_t< !std::is_void< _ud_ops_t>::value >
//   resetContainersToZero(){
//     for (std::size_t i=0; i<n; i++)
//       udOps_.get().set_zero(*data_(i).data());
//   }

// private:
//   typename std::conditional<
//   std::is_void<ud_ops_t>::value,
//   ::pressio::utils::impl::empty,
//   std::reference_wrapper<const ud_ops_t>
//   >::type udOps_;

//   std::reference_wrapper<const reconstuctor_type> fomStateReconstrObj_;

//   // data[0] contains the current fom state, i.e. step = n
//   // data[1] contains fom state at step n-1
//   // data[2] contains fom state at n-2
//   // etc..
//   data_type data_;
// };
#endif  // ROM_FOM_STATES_MANAGEMENT_ROM_MANAGER_FOM_STATES_STATIC_HPP_
