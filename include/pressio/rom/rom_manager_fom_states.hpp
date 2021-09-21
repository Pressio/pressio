/*
//@HEADER
// ************************************************************************
//
// rom_manager_fom_states.hpp
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

// =================================
//
// STEADY
//
// a notion of time stencil does not apply
// and we only need a single fom state
//
// =================================
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

  template <class RomStateType>
  void reconstructCurrentFomState(const RomStateType & romStateIn){
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


// =================================
//
// STATIC manager with stencil notion
//
// when we know the # of states at compile time
// here we currently only support the case
// where we have y_n+1, y_n, y_n-1, etc
//
// =================================
template <class FomStateType, class ReconstuctorType, std::size_t N>
class ManagerStencilFomStatesStatic
{
  static_assert(N>=1, "ManagerFomStates cannot be empty.");

public:
  using data_type  = std::array<FomStateType, N>;
  using value_type = FomStateType;

  ManagerStencilFomStatesStatic() = delete;
  ManagerStencilFomStatesStatic(const ManagerStencilFomStatesStatic &) = default;
  ManagerStencilFomStatesStatic & operator=(const ManagerStencilFomStatesStatic &) = delete;
  ManagerStencilFomStatesStatic(ManagerStencilFomStatesStatic &&) = default;
  ManagerStencilFomStatesStatic & operator=(ManagerStencilFomStatesStatic &&) = delete;
  ~ManagerStencilFomStatesStatic() = default;

  // constructor for n == 1
  template <std::size_t _N = N, mpl::enable_if_t<_N == 1, int> = 0>
  ManagerStencilFomStatesStatic(const ReconstuctorType & fomStateReconstr,
				const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 2
  template <std::size_t _N = N, mpl::enable_if_t<_N == 2, int> = 0>
  ManagerStencilFomStatesStatic(const ReconstuctorType & fomStateReconstr,
				const FomStateType & fomState)
    : fomStateReconstrObj_(fomStateReconstr),
      data_{::pressio::ops::clone(fomState),
            ::pressio::ops::clone(fomState)}
      {
        this->setZero();
      }

  // constructor for n == 3
  template <std::size_t _N = N, mpl::enable_if_t<_N == 3, int> = 0>
  ManagerStencilFomStatesStatic(const ReconstuctorType & fomStateReconstr,
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

  // n+1
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=1, FomStateType const &>
  operator()(::pressio::ode::nPlusOne) const {return data_[0];}

  // n
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=2, FomStateType const &>
  operator()(::pressio::ode::n) const {return data_[1];}

  // n-1
  template <std::size_t _N = N>
  mpl::enable_if_t<_N>=3, FomStateType const &>
  operator()(::pressio::ode::nMinusOne) const {return data_[2];}

  // ** methods to reconstruct fom state **
  template <class RomStateType>
  void reconstructCurrentFomState(const RomStateType & romStateIn)
  {
    this->reconstructAt(romStateIn, ::pressio::ode::nPlusOne());
  }

  // n+1
  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructAt(const RomStateType & romStateIn,
		::pressio::ode::nPlusOne)
  {
    fomStateReconstrObj_(romStateIn, data_[0]);
  }

  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N>=1 >
  reconstructFomStateAt(const RomStateType & romStateIn,
			::pressio::ode::nPlusOne tag)
  {
    reconstructAt(romStateIn, tag);
  }

  // n
  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
  reconstructAt(const RomStateType & romStateIn,
		::pressio::ode::n)
  {
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N>=2 >
  reconstructFomStateAt(const RomStateType & romStateIn,
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
  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N==2 >
  reconstructAtAndUpdatePrevious(const RomStateType & romStateIn,
				 ::pressio::ode::n /*tag*/)
  {
    /* when n == 2, it means I only have n+1 and n
     * so to reconstruct at n, I can simply
     * overwrite the data in data_[1] */
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  /* when n == 3, we have y_n+1, y_n, y_n-1 */
  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N==3>
  reconstructAtAndUpdatePrevious(const RomStateType & romStateIn,
				 ::pressio::ode::n /*tag*/)
  {
    /*
     * copy y_n into y_n-1
     * then reconstruct y_n
     */
    ::pressio::ops::deep_copy(data_[2], data_[1]);
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  /* when n == 4, we have y_n+1, y_n, y_n-1, y_n-2 */
  template <class RomStateType, std::size_t _N = N>
  mpl::enable_if_t< _N==4>
  reconstructAtAndUpdatePrevious(const RomStateType & romStateIn,
				 ::pressio::ode::n /*tag*/)
  {
    /*
     * copy y_n-1 into y_n-2
     * copy y_n into y_n-1
     * then reconstruct y_n */
    ::pressio::ops::deep_copy(data_[3], data_[2]);
    ::pressio::ops::deep_copy(data_[2], data_[1]);
    fomStateReconstrObj_(romStateIn, data_[1]);
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




// =================================
//
// DYNAMIC manager with stencil notion
//
// =================================
template <class FomStateType, class ReconstuctorType, class StencilEndsAt>
class ManagerStencilFomStatesDynamic;


//
// partial specialize for: n,n-1, n-2, etc
//
template <class FomStateType, class ReconstuctorType>
class ManagerStencilFomStatesDynamic<
  FomStateType, ReconstuctorType, ::pressio::ode::n>
{

public:
  using data_type  = std::vector<FomStateType>;
  using value_type = FomStateType;

  ManagerStencilFomStatesDynamic() = delete;
  ManagerStencilFomStatesDynamic(const ManagerStencilFomStatesDynamic &) = default;
  ManagerStencilFomStatesDynamic & operator=(const ManagerStencilFomStatesDynamic &) = delete;
  ManagerStencilFomStatesDynamic(ManagerStencilFomStatesDynamic &&) = default;
  ManagerStencilFomStatesDynamic & operator=(ManagerStencilFomStatesDynamic &&) = delete;
  ~ManagerStencilFomStatesDynamic() = default;

  ManagerStencilFomStatesDynamic(const ReconstuctorType & fomStateReconstr,
				   std::initializer_list<FomStateType> il)
    : fomStateReconstrObj_(fomStateReconstr), data_(il)
  {
    this->setZero();
  }

public:
  const std::size_t size() const { return data_.size(); }

  // n
  FomStateType const & operator()(::pressio::ode::n) const {
    assert(data_.size() >=1); return data_[0];
  }

  // n-1
  FomStateType const & operator()(::pressio::ode::nMinusOne) const {
    assert(data_.size() >=2);
    return data_[1];
  }

  template <class RomStateType>
  void reconstructCurrentFomState(const RomStateType & romStateIn)
  {
    this->reconstructAt(romStateIn, ::pressio::ode::n());
  }

  // n
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::n /*tag*/){
    assert(data_.size() >=1);
    fomStateReconstrObj_(romStateIn, data_[0]);
  }

  // n-1
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::nMinusOne /*tag*/){
    assert(data_.size() >=2);
    fomStateReconstrObj_(romStateIn, data_[1]);
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


//
// partial specialize for: n+1, n, n-1, n-2, etc
//
template <class FomStateType, class ReconstuctorType>
class ManagerStencilFomStatesDynamic<
  FomStateType, ReconstuctorType, ::pressio::ode::nPlusOne>
{

public:
  using data_type  = std::vector<FomStateType>;
  using value_type = FomStateType;

  ManagerStencilFomStatesDynamic() = delete;
  ManagerStencilFomStatesDynamic(const ManagerStencilFomStatesDynamic &) = default;
  ManagerStencilFomStatesDynamic & operator=(const ManagerStencilFomStatesDynamic &) = delete;
  ManagerStencilFomStatesDynamic(ManagerStencilFomStatesDynamic &&) = default;
  ManagerStencilFomStatesDynamic & operator=(ManagerStencilFomStatesDynamic &&) = delete;
  ~ManagerStencilFomStatesDynamic() = default;

  ManagerStencilFomStatesDynamic(const ReconstuctorType & fomStateReconstr,
				   std::initializer_list<FomStateType> il)
    : fomStateReconstrObj_(fomStateReconstr), data_(il)
  {
    this->setZero();
  }

public:
  const std::size_t size() const { return data_.size(); }

  // n+1
  FomStateType const & operator()(::pressio::ode::nPlusOne) const {
    assert(data_.size() >=1); return data_[0];
  }

  // n
  FomStateType const & operator()(::pressio::ode::n) const {
    assert(data_.size() >=2); return data_[1];
  }

  // n-1
  FomStateType const & operator()(::pressio::ode::nMinusOne) const {
    assert(data_.size() >=3); return data_[2];
  }

  // n-2
  FomStateType const & operator()(::pressio::ode::nMinusTwo) const {
    assert(data_.size() >=4); return data_[3];
  }

  // ** methods to reconstruct fom state **
  template <class RomStateType>
  void reconstructCurrentFomState(const RomStateType & romStateIn)
  {
    this->reconstructAt(romStateIn, ::pressio::ode::nPlusOne());
  }

  // n+1
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::nPlusOne /*tag*/){
    assert(data_.size() >=1);
    fomStateReconstrObj_(romStateIn, data_[0]);
  }

  // n
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::n /*tag*/){
    assert(data_.size() >=2);
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  // n-1
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::nMinusOne /*tag*/){
    assert(data_.size() >=3);
    fomStateReconstrObj_(romStateIn, data_[2]);
  }

  // n-2
  template <class RomStateType>
  void reconstructAt(const RomStateType & romStateIn, ::pressio::ode::nMinusTwo /*tag*/){
    assert(data_.size() >=4);
    fomStateReconstrObj_(romStateIn, data_[4]);
  }


  template <class RomStateType>
  void reconstructAtAndUpdatePrevious(const RomStateType & romStateIn,
				      ::pressio::ode::n /*tag*/)
  {

    assert(data_.size() >=2);

    if (data_.size() == 2){
      /* when n == 2, it means I only have n+1 and n
       * so reconstructing at n, we just overwrite data_[1] */
      fomStateReconstrObj_(romStateIn, data_[1]);
    }
    else if (data_.size() == 3){
      /* when n == 3, we have y_n+1, y_n, y_n-1 */
      // y_n becomes  y_n-1
      ::pressio::ops::deep_copy(data_[2], data_[1]);
      // reconstruct y_n
      fomStateReconstrObj_(romStateIn, data_[1]);
    }
    else if (data_.size() == 4){
      /* when n == 4, we have y_n+1, y_n, y_n-1, y_n-2 */

      // y_n-1 becomes y_n-2
      ::pressio::ops::deep_copy(data_[3], data_[2]);

      // y_n   becomes y_n-1
      ::pressio::ops::deep_copy(data_[2], data_[1]);

      // reconstruct y_n
      fomStateReconstrObj_(romStateIn, data_[1]);
    }
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
