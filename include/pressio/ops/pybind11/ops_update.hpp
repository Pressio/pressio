/*
//@HEADER
// ************************************************************************
//
// ops_update.hpp
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

#ifndef OPS_PYBIND11_OPS_RANK1_UPDATE_HPP_
#define OPS_PYBIND11_OPS_RANK1_UPDATE_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// V = a * V + b * V1
//----------------------------------------------------------------------
template<class T, class T1, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1>::value and
  mpl::variadic::any_of<_pybind_is_rank1, T, T1>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b)
{

  assert(extent(v,0) == extent(v1,0));
  assert(v.ndim() == v1.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type i=0; i<extent(v,0); ++i){
    v(i) = a*v(i) + b*v1(i);
  }
}

template<class T, class T1, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1>::value and
  mpl::variadic::any_of<_pybind_is_rank2, T, T1>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b)
{

  assert(v.ndim() == v1.ndim());
  assert(extent(v,0) == extent(v1,0));

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type j=0; j<extent(v,1); ++j){
    for (size_type i=0; i<extent(v,0); ++i){
      v(i,j) = a*v(i,j) + b*v1(i,j);
    }
  }
}

template<class T, class T1, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1>::value and
  mpl::variadic::all_of<_pybind_is_rank_dyn, T, T1>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b)
{

  assert(v.ndim() == v1.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  if (v.ndim()==1){
    assert(extent(v,0) == extent(v1,0));
    for (size_type i=0; i<extent(v,0); ++i){
      v(i) = a*v(i) + b*v1(i);
    }
  }
  else if (v.ndim()==2)
  {
    assert(extent(v,0) == extent(v1,0));
    assert(extent(v,1) == extent(v1,1));

    for (size_type j=0; j<extent(v,1); ++j){
      for (size_type i=0; i<extent(v,0); ++i){
	v(i,j) = a*v(i,j) + b*v1(i,j);
      }
    }
  }
}

//----------------------------------------------------------------------
// v = a*v + b*v1 + c*v2
//----------------------------------------------------------------------
template<class T, class T1, class T2, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2>::value and
  mpl::variadic::any_of<_pybind_is_rank1, T, T1, T2>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type i=0; i<extent(v,0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i);
  }
}

template<class T, class T1, class T2, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2>::value and
  mpl::variadic::any_of<_pybind_is_rank2, T, T1, T2>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(extent(v,1) == extent(v1,1));
  assert(extent(v,1) == extent(v2,1));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type j=0; j<extent(v,1); ++j){
    for (size_type i=0; i<extent(v,0); ++i){
      v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j);
    }
  }
}

template<class T, class T1, class T2, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2>::value and
  mpl::variadic::all_of<_pybind_is_rank_dyn, T, T1, T2>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c)
{

  using size_type = typename ::pressio::Traits<T>::size_type;
  if (v.ndim()==1){
    assert(extent(v,0) == extent(v1,0));
    assert(extent(v,0) == extent(v2,0));
    assert(v.ndim() == v1.ndim());
    assert(v.ndim() == v2.ndim());

    for (size_type i=0; i<extent(v,0); ++i){
      v(i) = a*v(i) + b*v1(i) + c*v2(i);
    }
  }
  else if (v.ndim()==2)
  {
    assert(extent(v,0) == extent(v1,0));
    assert(extent(v,0) == extent(v2,0));
    assert(extent(v,1) == extent(v1,1));
    assert(extent(v,1) == extent(v2,1));
    assert(v.ndim() == v1.ndim());
    assert(v.ndim() == v2.ndim());

    for (size_type j=0; j<extent(v,1); ++j){
      for (size_type i=0; i<extent(v,0); ++i){
	v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j);
      }
    }
  }
}


//----------------------------------------------------------------------
// v = a*v + b*v1 + c*v2 + d*v3
//----------------------------------------------------------------------
template<class T, class T1, class T2, class T3, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3>::value and
  mpl::variadic::any_of<_pybind_is_rank1, T, T1, T2, T3>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(extent(v,0) == extent(v3,0));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());
  assert(v.ndim() == v3.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type i=0; i<extent(v,0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i);
  }
}

template<class T, class T1, class T2, class T3, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3>::value and
  mpl::variadic::any_of<_pybind_is_rank2, T, T1, T2, T3>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(extent(v,0) == extent(v3,0));
  assert(extent(v,1) == extent(v1,1));
  assert(extent(v,1) == extent(v2,1));
  assert(extent(v,1) == extent(v3,1));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());
  assert(v.ndim() == v3.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type j=0; j<extent(v,1); ++j){
    for (size_type i=0; i<extent(v,0); ++i){
      v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j) + d*v3(i,j);
    }
  }
}

template<class T, class T1, class T2, class T3, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3>::value and
  mpl::variadic::all_of<_pybind_is_rank_dyn, T, T1, T2, T3>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d)
{

  assert(v.ndim() == v1.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  if (v.ndim()==1){
    assert(extent(v,0) == extent(v1,0));
    assert(extent(v,0) == extent(v2,0));
    assert(extent(v,0) == extent(v3,0));
    assert(v.ndim() == v1.ndim());
    assert(v.ndim() == v2.ndim());
    assert(v.ndim() == v3.ndim());

    for (size_type i=0; i<extent(v,0); ++i){
      v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i);
    }
  }
  else if (v.ndim()==2)
  {
    assert(extent(v,0) == extent(v1,0));
    assert(extent(v,0) == extent(v2,0));
    assert(extent(v,0) == extent(v3,0));
    assert(extent(v,1) == extent(v1,1));
    assert(extent(v,1) == extent(v2,1));
    assert(extent(v,1) == extent(v3,1));
    assert(v.ndim() == v1.ndim());
    assert(v.ndim() == v2.ndim());
    assert(v.ndim() == v3.ndim());

    for (size_type j=0; j<extent(v,1); ++j){
      for (size_type i=0; i<extent(v,0); ++i){
	v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j) + d*v3(i,j);
      }
    }
  }
}


//----------------------------------------------------------------------
// v = a*v + b*v1 + c*v2 + d*v3 + e*v4
//----------------------------------------------------------------------
template<class T, class T1, class T2, class T3, class T4, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3, T4>::value and
  mpl::variadic::any_of<_pybind_is_rank1, T, T1, T2, T3, T4>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d,
       const T4 & v4, scalar_t e)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(extent(v,0) == extent(v3,0));
  assert(extent(v,0) == extent(v4,0));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());
  assert(v.ndim() == v3.ndim());
  assert(v.ndim() == v4.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type i=0; i<extent(v,0); ++i){
    v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
  }
}

template<class T, class T1, class T2, class T3, class T4, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3, T4>::value and
  mpl::variadic::any_of<_pybind_is_rank2, T, T1, T2, T3, T4>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d,
       const T4 & v4, scalar_t e)
{

  assert(extent(v,0) == extent(v1,0));
  assert(extent(v,0) == extent(v2,0));
  assert(extent(v,0) == extent(v3,0));
  assert(extent(v,0) == extent(v4,0));
  assert(extent(v,1) == extent(v1,1));
  assert(extent(v,1) == extent(v2,1));
  assert(extent(v,1) == extent(v3,1));
  assert(extent(v,1) == extent(v4,1));
  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());
  assert(v.ndim() == v3.ndim());
  assert(v.ndim() == v4.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  for (size_type j=0; j<extent(v,1); ++j){
    for (size_type i=0; i<extent(v,0); ++i){
      v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j) + d*v3(i,j) + e*v4(i,j);
    }
  }
}

template<class T, class T1, class T2, class T3, class T4, class scalar_t>
::pressio::mpl::enable_if_t<
  mpl::variadic::all_of<_pybind_is_pybind, T, T1, T2, T3, T4>::value and
  mpl::variadic::all_of<_pybind_is_rank_dyn, T, T1, T2, T3, T4>::value
  >
update(T & v, scalar_t a,
       const T1 & v1, scalar_t b,
       const T2 & v2, scalar_t c,
       const T3 & v3, scalar_t d,
       const T4 & v4, scalar_t e)
{

  assert(v.ndim() == v1.ndim());
  assert(v.ndim() == v2.ndim());
  assert(v.ndim() == v3.ndim());
  assert(v.ndim() == v4.ndim());

  using size_type = typename ::pressio::Traits<T>::size_type;
  if (v.ndim()==1){
    for (size_type i=0; i<extent(v,0); ++i){
      v(i) = a*v(i) + b*v1(i) + c*v2(i) + d*v3(i) + e*v4(i);
    }
  }
  else if (v.ndim()==2)
  {
    for (size_type j=0; j<extent(v,1); ++j){
      for (size_type i=0; i<extent(v,0); ++i){
	v(i,j) = a*v(i,j) + b*v1(i,j) + c*v2(i,j) + d*v3(i,j) + e*v4(i,j);
      }
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_RANK1_UPDATE_HPP_
