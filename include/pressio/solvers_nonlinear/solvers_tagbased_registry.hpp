/*
//@HEADER
// ************************************************************************
//
// solvers_admissible_linear_solver_for_nonlinear_least_squares.hpp
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

#ifndef SOLVERS_NONLINEAR_TAG_BASED_REGISTRY_HPP_
#define SOLVERS_NONLINEAR_TAG_BASED_REGISTRY_HPP_


namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template<class, class> struct TagBasedStaticRegistry;

template<class ...Tags, class ...DataTypes>
struct TagBasedStaticRegistry< std::tuple<Tags...>, std::tuple<DataTypes...> >
{
  std::tuple<DataTypes...> d_;

  template<class ...CArgs>
  TagBasedStaticRegistry(CArgs && ... cargs) : d_(std::forward<CArgs>(cargs)...){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same, Tags...>::value)
      < mpl::size<Tags...>::value;
  }

  template<class TagToFind, class T>
  void set(T && o){
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    std::get<i>(d_) = std::forward<T>(o);
  }

  template<class TagToFind>
  auto & get(){
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    return std::get<i>(d_);
  }

  template<class TagToFind>
  const auto & get() const {
    constexpr int i = mpl::variadic::find_if_binary_pred_t<
      TagToFind, std::is_same, Tags...>::value;
    return std::get<i>(d_);
  }
};

template<class, class, class> struct TagBasedStaticRegistryExtension;

template<class Extendable, class ...ETags, class ...EDataTypes>
struct TagBasedStaticRegistryExtension<
  Extendable, std::tuple<ETags...>, std::tuple<EDataTypes...>
  >
{
  Extendable & reg_;
  using extension_registry_type = TagBasedStaticRegistry<
    std::tuple<ETags...>, std::tuple<EDataTypes...> >;
  extension_registry_type newReg_;

  template<class ...CArgs>
  TagBasedStaticRegistryExtension(Extendable & reg, CArgs && ... cargs)
    : reg_(reg), newReg_(std::forward<CArgs>(cargs)...){}

  template<class TagToFind>
  static constexpr bool contains(){
    return Extendable::template contains<TagToFind>() ||
      extension_registry_type::template contains<TagToFind>();
  }

  template<class TagToFind, class T>
  void set(T && o){
    if constexpr(Extendable::template contains<TagToFind>()){
      reg_.template set<TagToFind>(std::forward<T>(o));
    }
    else{
      newReg_.template set<TagToFind>(std::forward<T>(o));
    }
  }

  template<class TagToFind>
  auto & get(){
    if constexpr(Extendable::template contains<TagToFind>()){
      return reg_.template get<TagToFind>();
    }
    else{
      return newReg_.template get<TagToFind>();
    }
  }

  template<class TagToFind>
  const auto & get() const {
    if constexpr(Extendable::template contains<TagToFind>()){
      return reg_.template get<TagToFind>();
    }
    else{
      return newReg_.template get<TagToFind>();
    }
  }
};

template<
  class Tag, class DataType,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<Tag>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<Tag>, std::tuple<DataType>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2,
  class D1, class D2,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2>, std::tuple<D1,D2>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2, class T3,
  class D1, class D2, class D3,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T3>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2,T3>, std::tuple<D1,D2,D3>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

template<
  class T1, class T2, class T3, class T4,
  class D1, class D2, class D3, class D4,
  class Extendable, class ...CArgs
  >
auto reference_capture_registry_and_extend_with(Extendable & reg, CArgs && ... cargs)
{
  static_assert(!Extendable::template contains<T1>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T2>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T3>(),
		"Registry not extendable: already contains tag");
  static_assert(!Extendable::template contains<T4>(),
		"Registry not extendable: already contains tag");

  using ret_t = TagBasedStaticRegistryExtension<
    Extendable, std::tuple<T1,T2,T3,T4>, std::tuple<D1,D2,D3,D4>
    >;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

}}}
#endif  // SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_
