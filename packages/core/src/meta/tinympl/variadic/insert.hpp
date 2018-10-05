/*
//@HEADER
// ************************************************************************
//
//                                insert.hpp                               
//                         whatever
//              Copyright (C) 2015 Sandia Corporation
// This file was adapted from its original form in the tinympl library.
// The original file bore the following copyright:
//   Copyright (C) 2013, Ennio Barbaro.
// See LEGAL.md for more information.
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/


#ifndef TINYMPL_VARIADIC_INSERT_HPP
#define TINYMPL_VARIADIC_INSERT_HPP

#include <cstddef>

namespace tinympl {
namespace variadic {

/**
 * \ingroup VarBasics
 * \class insert
 * \brief Produce an output sequence from a variadic template inserting a new element at a given position
 * \param pos The position of the new element
 * \param T the new element
 * \param Out the output sequence type
 * \param Args... the input variadic template
 */
template<std::size_t pos,class T,template<class ... > class Out,class ... Args> struct insert;

template<std::size_t pos,class T,template<class ... > class Out,class Head,class ... Args> struct insert<pos,T,Out,Head,Args...>
{
private:
	template<class ... CopiedElements>
	struct impl
	{
		typedef typename insert<pos-1,T,Out,Args...>::template impl<CopiedElements...,Head>::type type;
	};

	template<std::size_t,class,template<class ... > class,class ...> friend struct insert;

public:
	static_assert(pos <= sizeof ... (Args) + 1,"pos > sequence size!");

	typedef typename impl<>::type type;
};

template<class T,template<class ... > class Out,class Head,class ... Args> struct insert<0,T,Out,Head,Args...>
{
private:
	template<class ... CopiedElements>
	struct impl
	{
		typedef Out< CopiedElements ..., T, Head,Args ... > type;
	};

	template<std::size_t,class,template<class ... > class,class ...> friend struct insert;

public:
	typedef typename impl<>::type type;
};

template<class T,template<class ... > class Out> struct insert<0,T,Out>
{
private:
	template<class ... CopiedElements>
	struct impl
	{
		typedef Out< CopiedElements ..., T > type;
	};

	template<std::size_t,class,template<class ... > class,class ...> friend struct insert;

public:
	typedef typename impl<>::type type;
};

} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_INSERT_HPP
