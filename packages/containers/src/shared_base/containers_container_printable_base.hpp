/*
//@HEADER
// ************************************************************************
//
// containers_container_printable_base.hpp
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

#ifndef CONTAINERS_SHARED_BASE_CONTAINER_PRINTABLE_BASE_HPP_
#define CONTAINERS_SHARED_BASE_CONTAINER_PRINTABLE_BASE_HPP_

#include "../containers_ConfigDefs.hpp"

namespace pressio{ namespace containers{

/*
 * number of elements to prints defaults = -1 , indicating to print all
 *
 * the char c = 'd', 'f'
 * 'd' = default, prints as they are, so for vectors prints vertically
 * 'f' = flatten, for vectors prints horizontally,
 *		  for matrices it falttens it
 */

template<typename derived_type, typename ord_t>
class ContainerPrintable1DBase
  : private utils::details::CrtpBase<
	ContainerPrintable1DBase<derived_type, ord_t>>{

  using this_t = ContainerPrintable1DBase<derived_type, ord_t>;
public:

  template <typename stream_t = std::ostream>
  void print(stream_t & os = std::cout,
	     char c = 'd',
	     ord_t numElementsToPrint = -1) const{
    this->underlying().printImpl(os, c, numElementsToPrint);
  }

  void printStdCout(char c = 'd', ord_t numElementsToPrint = -1) const{
    this->print(std::cout, c, numElementsToPrint);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<this_t>;

  ContainerPrintable1DBase() = default;
  ~ContainerPrintable1DBase() = default;
};//end class



template<typename derived_type, typename ord1_t, typename ord2_t = ord1_t>
class ContainerPrintable2DBase
  : private utils::details::CrtpBase<
  ContainerPrintable2DBase<derived_type, ord1_t, ord2_t>>{

  using this_t = ContainerPrintable2DBase<derived_type, ord1_t, ord2_t>;
public:

  template <typename stream_t = std::ostream>
  void print(stream_t & os = std::cout,
	     char c = 'd',
	     ord1_t niToPrint = -1,
	     ord2_t njToPrint = -1) const{
    this->underlying().printImpl(os, c, niToPrint, njToPrint);
  }

  void printStdCout(char c = 'd',
		    ord1_t niToPrint = -1,
		    ord2_t njToPrint = -1) const{
    this->print(std::cout, c, niToPrint, njToPrint);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<this_t>;

  ContainerPrintable2DBase() = default;
  ~ContainerPrintable2DBase() = default;

};//end class

}}//end namespace pressio::containers
#endif
