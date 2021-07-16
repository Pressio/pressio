/*
//@HEADER
// ************************************************************************
//
// solvers_base_observer.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OBSERVERS_SOLVERS_BASE_OBSERVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_OBSERVERS_SOLVERS_BASE_OBSERVER_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

struct BaseObserver
{
  using apply_function_type = void (*)(BaseObserver*, int, const void*);
  apply_function_type applyFnc_;

  BaseObserver() = default;
  BaseObserver(BaseObserver const &) = default;
  BaseObserver & operator=(BaseObserver const &) = default;
  BaseObserver(BaseObserver &&) = default;
  BaseObserver & operator=(BaseObserver &&) = default;
  virtual ~BaseObserver() = default;

  template <class StateType>
  void operator()(int step, const StateType & state)
  {
    (*applyFnc_)(this, step, &state);
  }
};

template <class StateType, class FunctorType>
class Observer : public BaseObserver
{
public:
  using state_type = StateType;

private:
  pressio::utils::instance_or_reference_wrapper<FunctorType> F_;

public:
  Observer(FunctorType Fin) : F_(Fin){}

  Observer() = delete;
  Observer(Observer const &) = default;
  Observer & operator=(Observer const &) = default;
  Observer(Observer &&) = default;
  Observer & operator=(Observer &&) = default;
  ~Observer() = default;

  FunctorType get() const{
    return F_.get();
  }
};

template<typename T>
void applyObserver(BaseObserver* observer,
		   int step,
		   const void * state_as_void)
{
  using state_t = typename T::state_type;
  const auto* p = static_cast<T*>(observer);
  const auto* state = reinterpret_cast<const state_t*>(state_as_void);
  p->get()(step, *state);
}

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OBSERVERS_SOLVERS_BASE_OBSERVER_HPP_
