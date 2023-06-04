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

#if 0
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
#endif


#define GETMETHOD(N) \
  template<class Tag, mpl::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  auto & get(){ return d##N##_; } \
  template<class Tag, mpl::enable_if_t< std::is_same<Tag, Tag##N >::value, int> = 0> \
  const auto & get() const { return d##N##_; }


template<class...> class TagBasedStaticRegistry;

template<class Tag1, class T1>
class TagBasedStaticRegistry<Tag1, T1>{
  T1 d1_;

public:
  template<class CA1>
  explicit TagBasedStaticRegistry(CA1 && a1) : d1_(std::forward<CA1>(a1)){}

  template<class TagToFind>
  static constexpr bool contains(){ return std::is_same<TagToFind, Tag1>::value; }

  GETMETHOD(1)
};

template<class Tag1, class Tag2, class T1, class T2>
class TagBasedStaticRegistry<Tag1, Tag2, T1, T2>{
  T1 d1_; T2 d2_;

public:
  template<class CA1, class CA2>
  TagBasedStaticRegistry(CA1 && a1, CA2 && a2)
    : d1_(std::forward<CA1>(a1)),
      d2_(std::forward<CA2>(a2)){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag2>::value) < 2;
  }
  GETMETHOD(1)
  GETMETHOD(2)
};

template<
  class Tag1, class Tag2, class Tag3, class Tag4, class Tag5, class Tag6,
  class T1, class T2, class T3, class T4, class T5, class T6
  >
class TagBasedStaticRegistry<
  Tag1, Tag2, Tag3, Tag4, Tag5, Tag6,
  T1, T2, T3, T4, T5, T6
  >
{
  T1 d1_; T2 d2_; T3 d3_; T4 d4_; T5 d5_; T6 d6_;

public:
  template<class CA1, class CA2, class CA3, class CA4, class CA5, class CA6>
  TagBasedStaticRegistry(CA1 && a1, CA2 && a2, CA3 && a3,
			 CA4 && a4, CA5 && a5, CA6 && a6)
    : d1_(std::forward<CA1>(a1)), d2_(std::forward<CA2>(a2)),
      d3_(std::forward<CA3>(a3)), d4_(std::forward<CA4>(a4)),
      d5_(std::forward<CA5>(a5)), d6_(std::forward<CA6>(a6)){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag2, Tag3, Tag4, Tag5, Tag6>::value) < 6;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
};

template<
  class Tag1, class Tag2, class Tag3, class Tag4, class Tag5, class Tag6, class Tag7, class Tag8,
  class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8
  >
class TagBasedStaticRegistry<
  Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8,
  T1, T2, T3, T4, T5, T6, T7, T8
  >
{
  T1 d1_; T2 d2_; T3 d3_; T4 d4_; T5 d5_; T6 d6_; T7 d7_; T8 d8_;

public:
  template<class CA1, class CA2, class CA3, class CA4, class CA5, class CA6, class CA7, class CA8>
  TagBasedStaticRegistry(CA1 && a1, CA2 && a2, CA3 && a3, CA4 && a4,
			 CA5 && a5, CA6 && a6, CA7 && a7, CA8 && a8)
    : d1_(std::forward<CA1>(a1)), d2_(std::forward<CA2>(a2)),
      d3_(std::forward<CA3>(a3)), d4_(std::forward<CA4>(a4)),
      d5_(std::forward<CA5>(a5)), d6_(std::forward<CA6>(a6)),
      d7_(std::forward<CA7>(a7)), d8_(std::forward<CA8>(a8)){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	   Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8>::value) < 8;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
};


template<
  class Tag1, class Tag2, class Tag3, class Tag4, class Tag5,
  class Tag6, class Tag7, class Tag8, class Tag9, class Tag10,
  class T1, class T2, class T3, class T4, class T5,
  class T6, class T7, class T8, class T9, class T10
  >
class TagBasedStaticRegistry<
 Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10,
  T1, T2, T3, T4, T5, T6, T7, T8, T9, T10
  >
{
  T1 d1_; T2 d2_; T3 d3_; T4 d4_; T5 d5_;
  T6 d6_; T7 d7_; T8 d8_; T9 d9_; T10 d10_;

public:
  template<
  class CA1, class CA2, class CA3, class CA4, class CA5,
  class CA6, class CA7, class CA8, class CA9, class CA10
  >
  TagBasedStaticRegistry(CA1 && a1, CA2 && a2, CA3 && a3, CA4 && a4,
			 CA5 && a5, CA6 && a6, CA7 && a7, CA8 && a8,
			 CA9 && a9, CA10 && a10)
    : d1_(std::forward<CA1>(a1)), d2_(std::forward<CA2>(a2)),
      d3_(std::forward<CA3>(a3)), d4_(std::forward<CA4>(a4)),
      d5_(std::forward<CA5>(a5)), d6_(std::forward<CA6>(a6)),
      d7_(std::forward<CA7>(a7)), d8_(std::forward<CA8>(a8)),
      d9_(std::forward<CA9>(a9)), d10_(std::forward<CA10>(a10)){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10>::value) < 10;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
};

template<
  class Tag1, class Tag2, class Tag3, class Tag4, class Tag5,
  class Tag6, class Tag7, class Tag8, class Tag9, class Tag10, class Tag11,
  class T1, class T2, class T3, class T4, class T5,
  class T6, class T7, class T8, class T9, class T10, class T11
  >
class TagBasedStaticRegistry<
  Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10, Tag11,
  T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11
  >
{
  T1 d1_; T2 d2_; T3 d3_; T4 d4_; T5 d5_;
  T6 d6_; T7 d7_; T8 d8_; T9 d9_; T10 d10_; T11 d11_;

public:
  template<
  class CA1, class CA2, class CA3, class CA4, class CA5,
  class CA6, class CA7, class CA8, class CA9, class CA10, class CA11
  >
  TagBasedStaticRegistry(CA1 && a1, CA2 && a2, CA3 && a3, CA4 && a4,
			 CA5 && a5, CA6 && a6, CA7 && a7, CA8 && a8,
			 CA9 && a9, CA10 && a10, CA11 && a11)
    : d1_(std::forward<CA1>(a1)), d2_(std::forward<CA2>(a2)),
      d3_(std::forward<CA3>(a3)), d4_(std::forward<CA4>(a4)),
      d5_(std::forward<CA5>(a5)), d6_(std::forward<CA6>(a6)),
      d7_(std::forward<CA7>(a7)), d8_(std::forward<CA8>(a8)),
      d9_(std::forward<CA9>(a9)), d10_(std::forward<CA10>(a10)),
      d11_(std::forward<CA11>(a11)){}

  template<class TagToFind>
  static constexpr bool contains(){
    return (mpl::variadic::find_if_binary_pred_t<TagToFind, std::is_same,
	    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8, Tag9, Tag10, Tag11>::value) < 11;
  }

  GETMETHOD(1)
  GETMETHOD(2)
  GETMETHOD(3)
  GETMETHOD(4)
  GETMETHOD(5)
  GETMETHOD(6)
  GETMETHOD(7)
  GETMETHOD(8)
  GETMETHOD(9)
  GETMETHOD(10)
  GETMETHOD(11)
};


template<class, class> struct TagBasedStaticRegistryTramp;

template<class ...Tags, class ...DataTypes>
struct TagBasedStaticRegistryTramp< std::tuple<Tags...>, std::tuple<DataTypes...> >{
  using type = TagBasedStaticRegistry<Tags..., DataTypes...>;
};

template<class TagsTuple, class DataTypesTuple>
using TagBasedStaticRegistryTramp_t =
  typename TagBasedStaticRegistryTramp<TagsTuple, DataTypesTuple>::type;

template<class Extendable, class ...Rest>
struct TagBasedStaticRegistryExtension
{
  Extendable & reg_;

  using extension_registry_type = TagBasedStaticRegistry<Rest...>;
  extension_registry_type newReg_;

  template<class ...CArgs>
  TagBasedStaticRegistryExtension(Extendable & reg, CArgs && ... cargs)
    : reg_(reg), newReg_(std::forward<CArgs>(cargs)...){}

  template<class TagToFind>
  static constexpr bool contains(){
    return Extendable::template contains<TagToFind>() ||
      extension_registry_type::template contains<TagToFind>();
  }

  // template<class TagToFind, class T>
  // mpl::enable_if_t< Extendable::template contains<TagToFind>() >
  // set(T && o){
  //   reg_.template set<TagToFind>(std::forward<T>(o));
  // }
  // template<class TagToFind, class T>
  // mpl::enable_if_t< !Extendable::template contains<TagToFind>() >
  // set(T && o){
  //   newReg_.template set<TagToFind>(std::forward<T>(o));
  // }

  template<
    class TagToFind,
    mpl::enable_if_t< Extendable::template contains<TagToFind>(), int> * = nullptr
    >
  auto & get(){
    return reg_.template get<TagToFind>();
  }

  template<
    class TagToFind,
    mpl::enable_if_t< !Extendable::template contains<TagToFind>(), int> * = nullptr
    >
  auto & get(){
    return newReg_.template get<TagToFind>();
  }

  template<
    class TagToFind,
    mpl::enable_if_t< Extendable::template contains<TagToFind>(), int> * = nullptr
    >
  const auto & get() const {
    return reg_.template get<TagToFind>();
  }

  template<
    class TagToFind,
    mpl::enable_if_t< !Extendable::template contains<TagToFind>(), int> * = nullptr
    >
  const auto & get() const {
    return newReg_.template get<TagToFind>();
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

  using ret_t = TagBasedStaticRegistryExtension<Extendable, Tag, DataType>;
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

  using ret_t = TagBasedStaticRegistryExtension<Extendable,T1,T2,D1,D2>;
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

  using ret_t = TagBasedStaticRegistryExtension<Extendable, T1,T2,T3, D1,D2,D3>;
  return ret_t(reg, std::forward<CArgs>(cargs)...);
}

}}}
#endif  // SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_
