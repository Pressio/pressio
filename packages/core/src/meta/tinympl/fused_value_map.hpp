
/*
//@HEADER
// ************************************************************************
//
//                           fused_value_map.hpp                           
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


#ifndef TINYMPL_FUSED_VALUE_MAP_HPP
#define TINYMPL_FUSED_VALUE_MAP_HPP

#include <algorithm>
#include "algorithm_variadic.hpp"

namespace tinympl {

template<class KeyType,class ValueType,KeyType... Keys>
struct fused_value_map
{
	static_assert( variadic::is_unique< std::integral_constant<KeyType,Keys>...>::type::value,"There are duplicate keys");

	typedef ValueType value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	// Capacity.
	static constexpr size_type size() noexcept { return sizeof ... (Keys); }
	static constexpr size_type max_size() { return size(); }
	static constexpr bool empty() noexcept { return size() == 0; }

	void fill(const value_type& __u) { std::fill_n(begin(), size(), __u); }
	void swap(fused_value_map& __other) noexcept(noexcept(swap(std::declval<reference>(), std::declval<reference>())))
	{ std::swap_ranges(begin(), end(), __other.begin()); }

	iterator begin() noexcept
	{ return iterator(data()); }

	const_iterator begin() const noexcept
	{ return const_iterator(data()); }

	iterator end() noexcept
	{ return iterator(data() + size()); }

	const_iterator end() const noexcept
	{ return const_iterator(data() + size()); }

	reverse_iterator rbegin() noexcept
	{ return reverse_iterator(end()); }

	const_reverse_iterator rbegin() const noexcept
	{ return const_reverse_iterator(end()); }

	reverse_iterator rend() noexcept
	{ return reverse_iterator(begin()); }

	const_reverse_iterator rend() const noexcept
	{ return const_reverse_iterator(begin()); }

	const_iterator cbegin() const noexcept
	{ return const_iterator(data()); }

	const_iterator cend() const noexcept
	{ return const_iterator(data() + size()); }

	const_reverse_iterator crbegin() const noexcept
	{ return const_reverse_iterator(end()); }

	const_reverse_iterator crend() const noexcept
	{ return const_reverse_iterator(begin()); }

	template<KeyType k>
	reference at() {
		enum {value = variadic::find<
			std::integral_constant<KeyType,k>,
			std::integral_constant<KeyType,Keys>...>::type::value};
		static_assert( value < size(), "Key not present in map");
		return elems_[value];
	}

	template<KeyType k>
	const_reference at() const {
		enum {value = variadic::find<
			std::integral_constant<KeyType,k>,
			std::integral_constant<KeyType,Keys>...>::type::value};
		static_assert( value < size(), "Key not present in map");
		return elems_[value];
	}

	pointer data() noexcept { return std::addressof( elems_[0] ); }
	const_pointer data() const noexcept { return std::addressof( elems_[0] ); }

	value_type elems_[ size() > 0 ? size() : 1 ];

	inline friend bool operator==(const fused_value_map & lhs, const fused_value_map & rhs)
    { return std::equal(lhs.begin(), lhs.end(), rhs.begin()); }

    inline friend bool operator!=(const fused_value_map & lhs, const fused_value_map & rhs)
    { return !(lhs == rhs); }

    inline friend bool operator<(const fused_value_map & lhs, const fused_value_map & rhs)
	{ return std::lexicographical_compare(lhs.begin(), lhs.end(),rhs.begin(), rhs.end()); }
	inline friend bool operator>(const fused_value_map & lhs, const fused_value_map & rhs)
	{ return rhs < lhs; }

	inline friend bool operator<=(const fused_value_map & lhs,const fused_value_map & rhs)
	{ return !(lhs > rhs); }

	inline friend bool operator>=(const fused_value_map & lhs,const fused_value_map & rhs)
	{ return !(lhs < rhs); }
};

}

#endif // TINYMPL_FUSED_VALUE_MAP_HPP
