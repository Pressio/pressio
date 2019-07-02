
#ifndef CONTAINERS_VECTOR_VECTOR_SHAREDMEM_BINARY_EXPRESSION_TEMPLATES_HPP_
#define CONTAINERS_VECTOR_VECTOR_SHAREDMEM_BINARY_EXPRESSION_TEMPLATES_HPP_

#include "../containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace exprtemplates{


template <typename der_t>
class SharedMemVecExpressionBase {
  SharedMemVecExpressionBase() = default;
  ~SharedMemVecExpressionBase() = default;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<der_t>::type;

public:
  static constexpr bool is_shared_mem = true;
  static constexpr bool is_vector_expression = true;
};
///-----------------------------------------------------

template <typename T,
	  typename enable = void>
struct is_sharedmem_vector_expression : std::false_type{};

template <typename T>
struct is_sharedmem_vector_expression<T,
      ::pressio::mpl::enable_if_t<
       ::pressio::mpl::publicly_inherits_from<
	T,SharedMemVecExpressionBase<T>
	 >::value
	>> : std::true_type{};


//----------------------------------------------------
// default
template <typename OP_t, typename T1,
	  typename T2, typename value_t,
	  typename ord_t, typename enable = void>
class SharedMemVectorBinaryExp
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>>{

  const T1 & a_ = {};
  const T2 & b_ = {};
  OP_t op_ = {};
  using this_t = SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>;
  friend SharedMemVecExpressionBase<this_t>;

 public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(T1 const& a, T2 const& b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_(i));}

  ord_t size() const {
    return a_.size();}
};


//----------------------------------------------------
// T1: whatever, T2: vector
template <typename OP_t, typename T1, typename T2,
	  typename value_t, typename ord_t>
class SharedMemVectorBinaryExp<
         OP_t, T1, T2, value_t, ord_t,
	 ::pressio::mpl::enable_if_t<
	   !std::is_scalar<T1>::value &&
	   containers::meta::is_vector_wrapper<T2>::value &&
	   containers::details::traits<T2>::is_shared_mem>
  >
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>>{

  const T1 & a_ = {};
  const T2 & b_ = {};
  OP_t op_ = {};
  using this_t = SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>;
  friend SharedMemVecExpressionBase<this_t>;

public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(T1 const& a, T2 const& b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_(i));}

  ord_t size() const {
    return b_.size();}
};


//-----------------------------------------------------
// T1: not scalar, T2: scalar
template <typename OP_t, typename T1,
	    typename value_t, typename ord_t>
class SharedMemVectorBinaryExp<
         OP_t, T1, value_t, value_t, ord_t,
	 ::pressio::mpl::enable_if_t<
	   !std::is_scalar<T1>::value &&
	   std::is_scalar<value_t>::value
	   > >
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,value_t,value_t,ord_t>>{

  const T1 & a_ = {};
  value_t b_ = {};
  OP_t op_ = {};
  using th_t = SharedMemVectorBinaryExp<OP_t,T1,value_t,value_t,ord_t>;
  friend SharedMemVecExpressionBase<th_t>;

public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(const T1 & a, value_t b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_);}

  ord_t size() const{
    return a_.size();}
};


}}}//end namespace pressio::containers::exprtemplates
#endif
