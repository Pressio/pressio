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

template <class FomStateType, class ReconstuctorType>
class ManagerFomStatesSteady
{
public:
  using data_type = std::array<FomStateType, 1>;
  using value_type = FomStateType;

  ManagerFomStatesSteady() = delete;
  ManagerFomStatesSteady(const ManagerFomStatesSteady &) = default;
  ManagerFomStatesSteady & operator=(const ManagerFomStatesSteady &) = delete;
  ManagerFomStatesSteady(ManagerFomStatesSteady &&) = default;
  ManagerFomStatesSteady & operator=(ManagerFomStatesSteady &&) = delete;
  ~ManagerFomStatesSteady() = default;

  ManagerFomStatesSteady(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState)}
  {
    this->setZero();
  }

public:
  static constexpr std::size_t size(){ return data_type::size(); }

  FomStateType const & currentFomState() const {return data_[0];}

  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romStateIn)
  {
    fomStateReconstrObj_(romStateIn, data_[0]);
  }

private:
  void setZero(){
    ::pressio::ops::set_zero(data_[0]);
  }

private:
  std::reference_wrapper<const ReconstuctorType> fomStateReconstrObj_;
  data_type data_;
};


// explicit stepping
template <class FomStateType, class ReconstuctorType, std::size_t N>
class ManagerFomStatesUnsteadyExplicit
{
  //
  // for ROMs with explicit time stepping we define indexing:
  // "n, n-1, n-2", etc
  //

  static_assert(N>=1, "ManagerFomStatesUnsteadyExplicit cannot be empty.");

public:
  using data_type  = std::array<FomStateType, N>;
  using value_type = FomStateType;

  ManagerFomStatesUnsteadyExplicit() = delete;
  ManagerFomStatesUnsteadyExplicit(const ManagerFomStatesUnsteadyExplicit &) = default;
  ManagerFomStatesUnsteadyExplicit & operator=(const ManagerFomStatesUnsteadyExplicit &) = delete;
  ManagerFomStatesUnsteadyExplicit(ManagerFomStatesUnsteadyExplicit &&) = default;
  ManagerFomStatesUnsteadyExplicit & operator=(ManagerFomStatesUnsteadyExplicit &&) = delete;
  ~ManagerFomStatesUnsteadyExplicit() = default;

  // constructor for n == 1
  template <std::size_t _N = N, mpl::enable_if_t<_N == 1, int> = 0>
  ManagerFomStatesUnsteadyExplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 2
  template <std::size_t _N = N, mpl::enable_if_t<_N == 2, int> = 0>
  ManagerFomStatesUnsteadyExplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 3
  template <std::size_t _N = N, mpl::enable_if_t<_N == 3, int> = 0>
  ManagerFomStatesUnsteadyExplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }


public:
  static constexpr std::size_t size(){ return data_type::size(); }

  // n
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=1, FomStateType const &>
  fomStateAt(::pressio::ode::n) const {return data_[0];}

  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=1, FomStateType const &>
  operator()(::pressio::ode::n) const {return data_[0];}

  // n-1
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=2, FomStateType const &>
  fomStateAt(::pressio::ode::nMinusOne) const {return data_[1];}

  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=2, FomStateType const &>
  operator()(::pressio::ode::nMinusOne) const {return data_[1];}

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
  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructAt(const rom_state_t & romStateIn, ::pressio::ode::n)
  {
    fomStateReconstrObj_(romStateIn, data_[0]);
  }

  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructFomStateAt(const rom_state_t & romStateIn,	::pressio::ode::n tag)
  {
    reconstructAt(romStateIn, tag);
  }

  // n-1
  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
  reconstructAt(const rom_state_t & romStateIn, ::pressio::ode::n)
  {
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
  reconstructFomStateAt(const rom_state_t & romStateIn, ::pressio::ode::n tag)
  {
    reconstructAt(romStateIn, tag);
  }

private:
  void setZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      ::pressio::ops::set_zero(data_[i]);
  }

private:
  std::reference_wrapper<const ReconstuctorType> fomStateReconstrObj_;
  data_type data_;
};


// for implicit stepping
template <class FomStateType, class ReconstuctorType, std::size_t N>
class ManagerFomStatesUnsteadyImplicit
{
  static_assert(N>=1, "ManagerFomStates cannot be empty.");

public:
  using data_type  = std::array<FomStateType, N>;
  using value_type = FomStateType;

  ManagerFomStatesUnsteadyImplicit() = delete;
  ManagerFomStatesUnsteadyImplicit(const ManagerFomStatesUnsteadyImplicit &) = default;
  ManagerFomStatesUnsteadyImplicit & operator=(const ManagerFomStatesUnsteadyImplicit &) = delete;
  ManagerFomStatesUnsteadyImplicit(ManagerFomStatesUnsteadyImplicit &&) = default;
  ManagerFomStatesUnsteadyImplicit & operator=(ManagerFomStatesUnsteadyImplicit &&) = delete;
  ~ManagerFomStatesUnsteadyImplicit() = default;

  // constructor for n == 1
  template <std::size_t _N = N, mpl::enable_if_t<_N == 1, int> = 0>
  ManagerFomStatesUnsteadyImplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 2
  template <std::size_t _N = N, mpl::enable_if_t<_N == 2, int> = 0>
  ManagerFomStatesUnsteadyImplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 3
  template <std::size_t _N = N, mpl::enable_if_t<_N == 3, int> = 0>
  ManagerFomStatesUnsteadyImplicit(const ReconstuctorType & fomStateReconstr,
                   const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

public:
  static constexpr std::size_t size(){ return data_type::size(); }

  // for implicit time stepping it makese sense to
  // index using "n+1, n, n-1, n-2", etc

  // ** methods to extract const ref to data **
  // n+1
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=1, FomStateType const &>
  fomStateAt(::pressio::ode::nPlusOne) const {return data_[0];}

  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=1, FomStateType const &>
  operator()(::pressio::ode::nPlusOne) const {return data_[0];}

  // n
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=2, FomStateType const &>
  fomStateAt(::pressio::ode::n) const {return data_[1];}

  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=2, FomStateType const &>
  operator()(::pressio::ode::n) const {return data_[1];}

  // n-1
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=3, FomStateType const &>
  fomStateAt(::pressio::ode::nMinusOne) const {return data_[2];}

  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=3, FomStateType const &>
  operator()(::pressio::ode::nMinusOne) const {return data_[2];}

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
  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::nPlusOne)
  {
    fomStateReconstrObj_.get()(romStateIn, data_[0]);
  }

  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructFomStateAt(const rom_state_t & romStateIn,
			::pressio::ode::nPlusOne tag)
  {
    reconstructAt(romStateIn, tag);
  }

  // n
  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
  reconstructAt(const rom_state_t & romStateIn,
		::pressio::ode::n)
  {
    fomStateReconstrObj_.get()(romStateIn, data_[1]);
  }

  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
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
  template <typename rom_state_t, std::size_t _N = N>
  mpl::enable_if_t< _N==2 >
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /* when n == 2, it means I only have n+1 and n
     * so to reconstruct at n, I can simply
     * overwrite the data in data_[1] */
    fomStateReconstrObj_.get()(romStateIn, data_[1]);
  }

  /* when n == 3, we have y_n+1, y_n, y_n-1 */
  template <
    typename rom_state_t,
    std::size_t _N = N
    >
  mpl::enable_if_t< _N==3>
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /*
     * copy y_n into y_n-1
     * then reconstruct y_n
     */
    ::pressio::ops::deep_copy(data_[2], data_[1]);
    fomStateReconstrObj_.get()(romStateIn, data_[1]);
  }

  /* when n == 4, we have y_n+1, y_n, y_n-1, y_n-2 */
  template <
    typename rom_state_t,
    std::size_t _N = N
    >
  mpl::enable_if_t< _N==4>
  reconstructWithStencilUpdate(const rom_state_t & romStateIn)
  {
    /*
     * copy y_n-1 into y_n-2
     * copy y_n into y_n-1
     * then reconstruct y_n */
    ::pressio::ops::deep_copy(data_[3], data_[2]);
    ::pressio::ops::deep_copy(data_[2], data_[1]);
    fomStateReconstrObj_.get()(romStateIn, data_[1]);
  }

private:
  void setZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      ::pressio::ops::set_zero(data_[i]);
  }

private:
  std::reference_wrapper<const ReconstuctorType> fomStateReconstrObj_;
  data_type data_;
};

}}//end namespace pressio::rom

#endif  // ROM_FOM_STATES_MANAGEMENT_ROM_MANAGER_FOM_STATES_STATIC_HPP_
