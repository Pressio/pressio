
/*
//@HEADER
// ************************************************************************
//
//                                ratio.hpp                                
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


#ifndef TINYMPL_RATIO_HPP
#define TINYMPL_RATIO_HPP

#include <cstdint>
#include <ratio>
#include "functional.hpp"

namespace tinympl
{

namespace detail
{
//Construct a rational number by forcing a coprime representation
template<std::intmax_t Num, std::intmax_t Den> struct make_rational
{
	typedef std::ratio<Num,Den> ratio_t;
	typedef std::ratio<
		ratio_t::num,
		ratio_t::den> type;
};
}

/**
 * \ingroup NewTypes
 * @{
 */

/**
 * \brief Convenience wrapper around `std::ratio` to automatically reduce `num` and `den` to coprime factors.
 *
 * `std::is_same< std::ratio<4,2>, std::ratio<2,1> >::value` is `false`, while
 * `std::is_same< rational<4,2>, rational<2,1> >::value` is `true`.
 *
 * `rational` forwards to `std::ratio`. The comparison functionals \ref plus, \ref minus, \ref multiplies, \ref divides,
 * \ref negate, \ref equal_to and \ref less are specialized to work transparently on `std::ratio`.
 */
template<std::intmax_t Num,std::intmax_t Den> using rational = typename detail::make_rational<Num,Den>::type;

template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct plus<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > : std::ratio_add<std::ratio<Num1,Den1>, std::ratio<Num2,Den2> > {};

template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct minus<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > : std::ratio_subtract<std::ratio<Num1,Den1>, std::ratio<Num2,Den2> > {};


template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct multiplies<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > : std::ratio_multiply<std::ratio<Num1,Den1>, std::ratio<Num2,Den2> > {};

template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct divides<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > : std::ratio_divide<std::ratio<Num1,Den1>, std::ratio<Num2,Den2> > {};

template<std::intmax_t Num,std::intmax_t Den> struct negate<std::ratio<Num,Den> >
{
	typedef std::ratio<-Num,Den> type;
};

template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct equal_to<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > :
	std::integral_constant<bool,
		std::ratio<Num1,Den1>::num == std::ratio<Num2,Den2>::num &&
		std::ratio<Num1,Den1>::den == std::ratio<Num2,Den2>::den> {};

template<std::intmax_t Num1,std::intmax_t Den1,
		std::intmax_t Num2,std::intmax_t Den2> struct less<
			std::ratio<Num1,Den1>,
			std::ratio<Num2,Den2> > : std::ratio_less<std::ratio<Num1,Den1>, std::ratio<Num2,Den2> > {};

/** @} */
}

#endif // TINYMPL_RATIO_HPP
