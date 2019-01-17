
#ifndef CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_
#define CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_SHAREDMEM_EIGEN_DYNAMIC_HPP_

#include "../base/core_multi_vector_sharedmem_base.hpp"

#include "../../meta/core_native_multi_vector_meta.hpp"
#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class MultiVector<wrapped_type,
		  meta::enable_if_t<
		    meta::is_dynamic_multi_vector_eigen<
		      wrapped_type>::value>
		  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorSharedMemBase< MultiVector<wrapped_type> >,
    public ContainerSubscriptable2DBase<
     MultiVector<wrapped_type>,
     typename details::traits<MultiVector<wrapped_type>>::scalar_t,
     typename details::traits<MultiVector<wrapped_type>>::ordinal_t>,
    public ContainerResizableBase<MultiVector<wrapped_type>, 2>
{

private:
  using this_t = MultiVector<wrapped_type>;
  using sc_t = typename details::traits<this_t>::scalar_t;
  using ord_t = typename details::traits<this_t>::ordinal_t;
  using wrap_t = typename details::traits<this_t>::wrapped_t;

public:
  MultiVector() = delete;

  explicit MultiVector(ord_t length, ord_t numVectors)
    : data_(length, numVectors){
        this->setZero();
  }

  explicit MultiVector(const wrap_t & other) : data_(other){}

  ~MultiVector() = default;

public:
  sc_t & operator()(ord_t irow, ord_t iVec){
    assert(iVec < this->numVectors() );
    assert(irow < this->length() );
    return data_(irow, iVec);
  }

  sc_t const & operator()(ord_t irow, ord_t iVec)const{
    assert(iVec < this->numVectors() );
    assert(irow < this->length() );
    return data_(irow, iVec);
  }

private:

  ord_t numVectorsImpl() const{
    return data_.cols();
  }

  ord_t lengthImpl() const {
    return data_.rows();
  };

  wrap_t * dataImpl(){
    return &data_;
  };

  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t dataCpImpl() const{
    return data_;
  };

  bool emptyImpl() const {
    return this->length()==0 ? true : false;
  }

  void scaleImpl(sc_t & factor){
    data_.coeffs() *= factor;
  };

  void setZeroImpl() {
    data_.setConstant(static_cast<sc_t>(0));
  }

  void resizeImpl(ord_t newlength, ord_t nVec){
    data_.resize(newlength, nVec);
    this->setZero();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorSharedMemBase< this_t >;
  friend ContainerSubscriptable2DBase< this_t, sc_t, ord_t>;
  friend ContainerResizableBase<this_t, 2>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::core
#endif
