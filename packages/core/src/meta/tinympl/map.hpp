
/*
//@HEADER
// ************************************************************************
//
//                                 map.hpp                                 
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


#ifndef TINYMPL_MAP_HPP
#define TINYMPL_MAP_HPP

#include "equal_to.hpp"
#include "variadic/all_of.hpp"
#include "variadic/is_unique.hpp"
#include "variadic/find.hpp"
#include "variadic/at.hpp"
#include "variadic/left_fold.hpp"
#include "variadic/remove_if.hpp"
#include "unordered_equal.hpp"
#include "sequence.hpp"
#include <type_traits>
#include <utility>

namespace tinympl {

/**
 * \brief Determine whether a type is an `std::pair`
 */
template<class T> struct is_pair : std::false_type {};

template<class FirstType,class SecondType>
struct is_pair<std::pair<FirstType,SecondType> > : std::true_type
{
	typedef FirstType first_type;
	typedef SecondType second_type;
};

/**
 * \ingroup Containers
 * @{
 */

/**
 * \class map
 * \brief A compile time map from a type to another
 * This class represents a compile time mapping between types. The mapping is specified using
 * `std::pair`, that is every parameter of `map` must be an `std::pair<KeyType,ValueType>`.
 *
 * `map` supports standard insertion/removal and element access.
 * \note \ref equal_to is specialized for `map`, in order to make the order in which the key/value pairs are specified irrelevant.
 */
template<class ... Args> struct map
{
	static_assert( variadic::all_of< is_pair, Args...>::type::value, "All the arguments of a map must be key/value pairs");
	static_assert( variadic::is_unique<typename Args::first_type ...>::type::value,"Duplicate keys in the map");

	/**
	 * \class at
	 * \brief Return the value element with the given key
	 */
	template<class T>
	struct at
	{
		enum {index = variadic::find<T, typename Args::first_type ... >::type::value };

		static_assert(index < sizeof...(Args),"Key T not present in the map");

		typedef typename variadic::at<index,Args...>::type::second_type type;
	};

	enum
	{
		size = sizeof ... (Args) //!< The number of elements contained in the map
	};

	enum
	{
		empty = (size == 0) //!< Determines whether the map is empty
	};

	/**
	 * \brief Count the number of elements in the map with a given key.
	 * \note Since this class is a `map` and not a `multimap`, the only possible results for this operation are 0 and 1.
	 */
	template<class Key>
	using count = std::integral_constant<
			std::size_t,
			(variadic::find<Key, typename Args::first_type ... >::type::value == size ? 0 : 1)>;

	/**
	 * \class insert
	 * \brief Returns a new map with the new Key/Value pair, or this map if the key is already present in the map
	 */
	template<class Key,class Value>
	struct insert : std::conditional<
			count<Key>::type::value == 0,
			map<Args..., std::pair<Key,Value> >,
			map<Args...> > {};

	/**
	 * \class insert_many
	 * \brief Calls \ref insert many times to insert many Key/Value pairs
	 */
	template<class ... KeyValuePairs>
	struct insert_many
	{
		static_assert( variadic::all_of< is_pair, KeyValuePairs...>::type::value, "All the arguments of insert_many must be key/value pairs");
		template<class Map,class T> using insert_one_t = typename Map::template insert<typename T::first_type, typename T::second_type>;

		typedef typename variadic::left_fold<insert_one_t,map,KeyValuePairs...>::type type;
	};

	/**
	 * \class erase
	 * \brief Return a new map constructed by the current map removing the `Key` key, if present, otherwise return the current map.
	 */
	template<class Key>
	class erase
	{
		template<class T> using key_comparer = std::is_same<typename T::first_type,Key>;

	public:
		typedef typename variadic::remove_if< key_comparer, tinympl::map, Args...>::type type;
	};
};

/** @} */

template<class ... As,class ... Bs>
struct equal_to< map<As...>, map<Bs...> > : unordered_equal< sequence<As...>, sequence<Bs...> > {};

}

#endif // TINYMPL_MAP_HPP
