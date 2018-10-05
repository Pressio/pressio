/*
//@HEADER
// ************************************************************************
//
//                          void_t.hpp
//                         tinympl
//              Copyright (C) 2017 NTESS, LLC
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

#ifndef TINYMPL_VOID_T_H_
#define TINYMPL_VOID_T_H_

namespace tinympl {

// A void_t implementation that works with gcc-4.9 (workaround for bug 64395
// From: http://stackoverflow.com/questions/35753920/why-does-the-void-t-detection-idiom-not-work-with-gcc-4-9
namespace _void_t_impl {

template <class... >
struct make_void { using type = void; };


} // end namepace _void_t_impl

template <class... T>
using void_t = typename _void_t_impl::make_void<T...>::type;

namespace _void_t_impl {

template <template <class...> class F>
struct make_void_template { using type = void; };

} // end namepace _void_t_impl

template <template <class...> class F>
using void_template_t = typename _void_t_impl::make_void_template<F>::type;


} // end namespace tinympl

#endif /* TINYMPL_VOID_T_H_ */
