
/*
//@HEADER
// ************************************************************************
//
//                              fused_map.hpp                              
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


#ifndef TINYMPL_FUSED_MAP_HPP
#define TINYMPL_FUSED_MAP_HPP

#include "variadic/all_of.hpp"
#include "variadic/is_unique.hpp"
#include "map.hpp"
#include <algorithm>

namespace tinympl {

template<class ... KeyValuePairs>
struct fused_map : std::tuple<typename KeyValuePairs::second_type ... >
{
	static_assert( variadic::all_of< is_pair, KeyValuePairs...>::type::value, "All the arguments of a map must be key/value pairs");
	static_assert( variadic::is_unique<typename KeyValuePairs::first_type ...>::type::value,"Duplicate keys in the map");

	typedef std::tuple<typename KeyValuePairs::second_type ... > base_type;
	typedef map<KeyValuePairs...> map_type;

	using base_type::base_type;

	template<class Key>
	typename map_type::template at<Key>::type & at()
	{
		return std::get< map_type::template at<Key>::index >(*this);
	}

	template<class Key>
	typename map_type::template at<Key>::type const & at() const
	{
		return std::get< map_type::template at<Key>::index >(*this);
	}

	enum {size = map_type::size};
	enum {empty = map_type::empty};

	template<class Key>
	using count = typename map_type::template count<Key>;
};

}

#endif // TINYMPL_FUSED_MAP_HPP
