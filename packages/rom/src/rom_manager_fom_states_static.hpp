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

#ifndef ROM_ROM_MANAGER_FOM_STATES_STATIC_HPP_
#define ROM_ROM_MANAGER_FOM_STATES_STATIC_HPP_

namespace pressio{ namespace rom{

template <typename fom_state_type, std::size_t n, typename reconstuctor_type, typename ud_ops_t>
class ManagerFomStatesStatic
{
  static_assert
  ( ::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
	"Currently, you can only create a ManagerFomStatesStatic of pressio wrappers.");

  static_assert
  (n>=1,
	"You are trying to instantiate a state FomStatesContainer object with zero size.\
   Something is not right, because this object would then be unusable.");

public:
  using data_type  = ::pressio::containers::IndexableStaticCollection<fom_state_type, n>;
  using value_type = fom_state_type;

  ManagerFomStatesStatic() = delete;
  ManagerFomStatesStatic(const ManagerFomStatesStatic &) = default;
  ManagerFomStatesStatic & operator=(const ManagerFomStatesStatic &) = default;
  ManagerFomStatesStatic(ManagerFomStatesStatic &&) = default;
  ManagerFomStatesStatic & operator=(ManagerFomStatesStatic &&) = default;
  ~ManagerFomStatesStatic() = default;

  template <typename ... Args>
  ManagerFomStatesStatic(const reconstuctor_type & fomStateReconstr,
			 Args && ... args)
    : fomStateReconstrObj_(fomStateReconstr),
      data_( std::forward<Args>(args)... )
  {
    this->resetContainersToZero();
  }

  template <typename ... Args>
  ManagerFomStatesStatic(const reconstuctor_type & fomStateReconstr,
			 const ud_ops_t * udOps,
			 Args && ... args)
    : udOps_(udOps),
      fomStateReconstrObj_(fomStateReconstr),
      data_( std::forward<Args>(args)... )
  {
    this->resetContainersToZero();
  }

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  const fom_state_type & getCRefToCurrentFomState() const{
    static_assert( n>=1, "Cannot call getCRefToCurrentFomState if n < 1");
    return data_(0);
  }

  const fom_state_type & getCRefToFomStatePrevStep() const{
    static_assert( n>=2, "Cannot call getCRefToFomStatePrevStep if n < 2");
    return data_(1);
  }

  const fom_state_type & getCRefToFomStatePrevPrevStep() const{
    static_assert( n>=3, "Cannot call getCRefToFomStatePrevPrevStep if n < 3");
    return data_(2);
  }


  /* this method reconstructs the current FOM state */
  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romStateIn)
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom state");
#endif

    static_assert( n>=1, "Cannot call reconstructCurrentFomState if n < 1");
    fomStateReconstrObj_.get()(romStateIn, data_(0));

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom state");
#endif
  }


  // overload the left shift operator to use when we need to
  // reconstruct the FOM state at n-1 and shift all previous ones

  /* when n == 1, disable the << operator since that is meant
   * to use for time-dep problems when we have to store the states
   * for the stepper stencil */

  /* when n == 2, it means I only have current state and previous one
   * so when I need to reconstruct the previous state, I can simply
   * overwrite the data in data_(1) */
  template < typename rom_state_t, std::size_t _n = n>
  ::pressio::mpl::enable_if_t<_n == 2>
  operator << (const rom_state_t & romStateIn)
  {
    // reconstrct the FOM state at n-1
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  /* when n >= 3, we need to deep copy data to
   * such that y_t-2 goes into y_t-3,
   * and y_t-1 goes into y_t-2, etc.
   * and then finally we overwrite data_(0) */
  template <typename _ud_ops_t = ud_ops_t, typename rom_state_t, std::size_t _n = n >
  ::pressio::mpl::enable_if_t< _n >= 3 and std::is_void<_ud_ops_t>::value >
  operator << (const rom_state_t & romStateIn)
  {
    // copy all states back, such that y_t-2 goes into y_t-3,
    // and y_t-1 goes into y_t-2, etc. so that y_t-1 is free to overwrite
    for (std::size_t i=n-2; i>=1; --i){
      const auto & src  = data_(i);
      auto & dest = data_(i+1);
      ::pressio::ops::deep_copy(dest, src);
    }
    // then, reconstrct the FOM state at t-1
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

  template <typename _ud_ops_t = ud_ops_t, typename rom_state_t, std::size_t _n = n >
    ::pressio::mpl::enable_if_t< _n >= 3 and !std::is_void<_ud_ops_t>::value >
  operator << (const rom_state_t & romStateIn)
  {
    for (std::size_t i=n-2; i>=1; --i){
      const auto & src  = data_(i);
      auto & dest = data_(i+1);
      udOps_->deep_copy(*dest.data(), *src.data());
    }
    // then, reconstrct the FOM state at t-1
    fomStateReconstrObj_.get()(romStateIn, data_(1));
  }

private:
  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void< _ud_ops_t>::value >
  resetContainersToZero(){
    for (std::size_t i=0; i<n; i++)
      ::pressio::ops::set_zero(data_(i));
  }

  template <typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void< _ud_ops_t>::value >
  resetContainersToZero(){
    for (std::size_t i=0; i<n; i++)
      udOps_->set_zero(*data_(i).data());
  }

private:
  const ud_ops_t * udOps_ = nullptr;
  std::reference_wrapper<const reconstuctor_type> fomStateReconstrObj_;

  // data[0] contains the current fom state, i.e. step = n
  // data[1] contains fom state at step n-1
  // data[2] contains fom state at n-2
  // etc..
  data_type data_;
};

}}//end namespace pressio::rom
#endif  // ROM_ROM_MANAGER_FOM_STATES_STATIC_HPP_
