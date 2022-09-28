
#ifndef PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
#define PRESSIO_ROM_LINEAR_SUBSPACE_HPP_

namespace pressio{ namespace rom{

template <class BasisMatrixType>
class LinearSubspace
{

public:
  enum class SpanningSet{Columns, Rows};
  using basis_matrix_type = std::remove_cv_t<BasisMatrixType>;

private:
  static_assert(std::is_class<basis_matrix_type>::value,
		"std::is_class<basis_matrix_type> == false");
  static_assert(std::is_copy_constructible< basis_matrix_type >::value,
		"template argument must be copy constructible");
  static_assert( !mpl::is_std_shared_ptr<basis_matrix_type>::value,
		"std::unique_ptr is not valid template arguments");
  static_assert(all_have_traits<basis_matrix_type>::value,
		"all_have_traits<basis_matrix_type>::value == false");
  static_assert(::pressio::Traits<basis_matrix_type>::rank == 2,
		"::pressio::Traits<basis_matrix_type>::rank != 2");

private:
  const basis_matrix_type basisMatrix_;
  const SpanningSet spanningSetValue_;

public:
  LinearSubspace(const basis_matrix_type & basisMatrix,
		 SpanningSet spanningSetValue)
    : basisMatrix_(::pressio::ops::clone(basisMatrix)),
      spanningSetValue_(spanningSetValue){}

  LinearSubspace(basis_matrix_type && basisMatrix,
		 SpanningSet spanningSetValue)
    : basisMatrix_(std::move(basisMatrix)),
      spanningSetValue_(spanningSetValue){}

  LinearSubspace(const LinearSubspace & other)
    : basisMatrix_(::pressio::ops::clone(other.basisMatrix_)),
      spanningSetValue_(other.spanningSetValue_){}

  /* For this class to be really immutable, we should not have a
     move assign operator and copy assign.
     Since we have const members, the compiler defines all those as deleted.
     And since we have a copy constructor, the move constructor does not
     particupare in OR so the copy constructor is always called.
  */

  ~LinearSubspace() = default;

  const basis_matrix_type & basis() const{
    return basisMatrix_;
  }

  const std::size_t dimension() const{
    switch(spanningSetValue_){
    case SpanningSet::Columns:
      return ::pressio::ops::extent(basisMatrix_, 1);
    case SpanningSet::Rows:
      return ::pressio::ops::extent(basisMatrix_, 0);
    default:
      return 0;
    }
  }

  bool isColumnSpace() const{
    return spanningSetValue_ == SpanningSet::Columns;
  }

  bool isRowSpace() const{
    return spanningSetValue_ == SpanningSet::Rows;
  }
};

}}
#endif
