/*
//@HEADER
// ************************************************************************
//
// containers_tensor_sharedmem_pybind.hpp
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

#ifndef CONTAINERS_TENSOR_CONCRETE_CONTAINERS_TENSOR_SHAREDMEM_PYBIND_HPP_
#define CONTAINERS_TENSOR_CONCRETE_CONTAINERS_TENSOR_SHAREDMEM_PYBIND_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Tensor<
  1, wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_array_pybind<wrapped_type>::value
    >
  >
{
public:
  using this_t	    = Tensor<1, wrapped_type>;
  using traits	    = details::traits<this_t>;
  using sc_t	    = typename traits::scalar_t;
  using ord_t	    = typename traits::ordinal_t;
  using wrap_t	    = typename traits::wrapped_t;
  using ref_t	    = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;
  using mut_proxy_t = typename traits::mut_proxy_t;
  using proxy_t	    = typename traits::proxy_t;
  static constexpr int rank = 1;

public:
  Tensor() = default;

  // constructor with extents
  template<int _rank = rank, mpl::enable_if_t<_rank==1, int> = 0>
  explicit Tensor(std::size_t ext) : data_(ext){}

  // constructor from pybind::array
  explicit Tensor(const wrap_t & src)
    : data_{ wrap_t(const_cast<wrap_t &>(src).request()) }
  {
    if ( data_.ndim() != 1 ){
      throw std::runtime_error
	("You cannot construct a rank-1 tensor wrapper from a pybind::array with ndim!=1");
    }
    const auto srcPx = src.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // -----------------------------
  // constructor with view semantics
  // use only this if you know what you are doing
  // -----------------------------
  Tensor(wrap_t src, ::pressio::view) : data_{src}
  {
    if ( data_.ndim() != 1 ){
      throw std::runtime_error
	("You cannot construct a rank-1 tensor wrapper from a pybind::array with ndim!=1");
    }
  }

  // copy constructor
  Tensor(Tensor const & other) : data_({ other.extent(0) })
  {
    const auto srcPx = other.data_.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // copy assignment
  Tensor & operator=(const Tensor & other) = delete;
  Tensor(Tensor && other) = default;
  Tensor & operator=(Tensor && o) = delete;
  ~Tensor(){};

public:
  ord_t extent(ord_t i) const { return data_.shape(i); }
  wrap_t const * data() const{ return &data_; }
  wrap_t * data(){ return &data_; }
  proxy_t proxy() const{ return data_.unchecked(); }
  mut_proxy_t proxy(){ return data_.mutable_unchecked(); }

  template<int _rank = rank>
  mpl::enable_if_t<_rank==1, ref_t>
  operator()(ord_t i){ return data_(i); };

  template<int _rank = rank>
  mpl::enable_if_t<_rank==1, const_ref_t>
  operator()(ord_t i) const{ return data_(i); };

private:
  template<typename source_proxy_t>
  void copySourceProxyToMine(source_proxy_t srcPx)
  {
    auto myProxy = data_.mutable_unchecked();
    for (ord_t i=0; i<srcPx.shape(0); ++i){
      myProxy(i) = srcPx(i);
    }
  }

private:
  wrap_t data_ = {};
};


template <typename wrapped_type>
class Tensor<
  2, wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_fstyle_array_pybind<wrapped_type>::value
    >
  >
{
public:
  using this_t	    = Tensor<2, wrapped_type>;
  using traits	    = details::traits<this_t>;
  using sc_t	    = typename traits::scalar_t;
  using ord_t	    = typename traits::ordinal_t;
  using wrap_t	    = typename traits::wrapped_t;
  using ref_t	    = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;
  using mut_proxy_t = typename traits::mut_proxy_t;
  using proxy_t	    = typename traits::proxy_t;
  static constexpr int rank = 2;

public:
  Tensor() = default;

  // constructor with extents
  template<int _rank = rank, mpl::enable_if_t<_rank==2, int> = 0>
  Tensor(std::size_t ext1, std::size_t ext2)
    : data_({ext1, ext2})
  {}

  // constructor from pybind::array with f-layout
  explicit Tensor(const wrap_t & src)
    : data_{ wrap_t(const_cast<wrap_t &>(src).request()) }
  {
    if ( data_.ndim() != 2 ){
      throw std::runtime_error
	("You cannot construct a rank-2 tensor wrapper from a pybind::array with ndim!=2");
    }
    const auto srcPx = src.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // -----------------------------
  // constructor with view semantics
  // use only this if you know what you are doing
  // -----------------------------
  Tensor(wrap_t src, ::pressio::view)
    : data_{src}
  {
    if ( data_.ndim() != 2 ){
      throw std::runtime_error
	("You cannot construct a rank-2 tensor wrapper from a pybind::array with ndim!=2");
    }
  }

  // copy constructor
  Tensor(Tensor const & other) : data_({ other.extent(0), other.extent(1) })
  {
    const auto srcPx = other.data_.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // copy assignment
  Tensor & operator=(const Tensor & other) = delete;
  Tensor(Tensor && other) = default;
  Tensor & operator=(Tensor && o) = delete;
  ~Tensor(){};

public:
  ord_t extent(ord_t i) const { return data_.shape(i); }
  wrap_t const * data() const{ return &data_; }
  wrap_t * data(){ return &data_; }
  proxy_t proxy() const{ return data_.unchecked(); }
  mut_proxy_t proxy(){ return data_.mutable_unchecked(); }

  template<int _rank = rank>
  mpl::enable_if_t<_rank==2, ref_t>
  operator()(ord_t i, ord_t j){ return data_(i,j); };

  template<int _rank = rank>
  mpl::enable_if_t<_rank==2, const_ref_t>
  operator()(ord_t i, ord_t j) const{ return data_(i,j); };

private:
  // these copy are best for column-major layout
  // should specialize for when wrapped_t is row major
  template<typename source_proxy_t>
  void copySourceProxyToMine(source_proxy_t srcPx)
  {
    auto myProxy = data_.mutable_unchecked();
    for (ord_t j=0; j<srcPx.shape(1); ++j){
      for (ord_t i=0; i<srcPx.shape(0); ++i){
	myProxy(i,j) = srcPx(i,j);
      }
    }
  }

private:
  wrap_t data_ = {};
};

template <typename wrapped_type>
class Tensor<
  3, wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_fstyle_array_pybind<wrapped_type>::value
    >
  >
{
public:
  using this_t	    = Tensor<3, wrapped_type>;
  using traits	    = details::traits<this_t>;
  using sc_t	    = typename traits::scalar_t;
  using ord_t	    = typename traits::ordinal_t;
  using wrap_t	    = typename traits::wrapped_t;
  using ref_t	    = typename traits::reference_t;
  using const_ref_t = typename traits::const_reference_t;
  using mut_proxy_t = typename traits::mut_proxy_t;
  using proxy_t	    = typename traits::proxy_t;
  static constexpr int rank = 3;

public:
  Tensor() = default;

  // constructor with extents
  template<int _rank = rank, mpl::enable_if_t<_rank==3, int> = 0>
  Tensor(std::size_t ext1, std::size_t ext2, std::size_t ext3)
    : data_({ext1, ext2, ext3})
  {}

  // constructor from pybind::array
  explicit Tensor(const wrap_t & src)
    : data_{ wrap_t(const_cast<wrap_t &>(src).request()) }
  {
    if ( data_.ndim() != 3 ){
      throw std::runtime_error
	("You cannot construct a rank-3 tensor wrapper from a pybind::array with ndim!=3");
    }

    const auto srcPx = src.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // -----------------------------
  // constructor with view semantics
  // use only this if you know what you are doing
  // -----------------------------
  Tensor(wrap_t src, ::pressio::view) : data_{src}
  {
    if ( data_.ndim() != 3 ){
      throw std::runtime_error
	("You cannot construct a rank-3 tensor wrapper from a pybind::array with ndim!=3");
    }
  }

  // copy constructor
  Tensor(Tensor const & other)
    : data_({ other.extent(0), other.extent(1), other.extent(2) })
  {
    const auto srcPx = other.data_.unchecked();
    copySourceProxyToMine(srcPx);
  }

  // copy assignment
  Tensor & operator=(const Tensor & other) = delete;
  Tensor(Tensor && other) = default;
  Tensor & operator=(Tensor && o) = delete;
  ~Tensor(){};

public:
  ord_t extent(ord_t i) const { return data_.shape(i); }
  wrap_t const * data() const{ return &data_; }
  wrap_t * data(){ return &data_; }
  proxy_t proxy() const{ return data_.unchecked(); }
  mut_proxy_t proxy(){ return data_.mutable_unchecked(); }

  template<int _rank = rank>
  mpl::enable_if_t<_rank==3, ref_t>
  operator()(ord_t i, ord_t j, ord_t k){ return data_(i,j,k); };

  template<int _rank = rank>
  mpl::enable_if_t<_rank==3, const_ref_t>
  operator()(ord_t i, ord_t j, ord_t k) const{ return data_(i,j, k); };

private:
  // these copy are best for column-major layout
  // should specialize for when wrapped_t is row major
  template<typename source_proxy_t>
  void copySourceProxyToMine(source_proxy_t srcPx)
  {
    auto myProxy = data_.mutable_unchecked();
    for (ord_t k=0; k<srcPx.shape(2); ++k){
      for (ord_t j=0; j<srcPx.shape(1); ++j){
	for (ord_t i=0; i<srcPx.shape(0); ++i){
	  myProxy(i,j,k) = srcPx(i,j,k);
	}
      }
    }
  }

private:
  wrap_t data_ = {};
};

}} //end namespace pressio::containers
#endif  // CONTAINERS_TENSOR_CONCRETE_CONTAINERS_TENSOR_SHAREDMEM_PYBIND_HPP_
