
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
  static_assert(std::is_copy_constructible< basis_matrix_type >::value,
		"template argument must be copy constructible");
  static_assert( !std::is_pointer<basis_matrix_type>::value,
		"basis_matrix_type cannot be pointers");
  static_assert( !mpl::is_std_shared_ptr<basis_matrix_type>::value,
		"std::unique_ptr is not valid template arguments");

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

  LinearSubspace& operator=(const LinearSubspace & /*other*/) = delete;

  /* For this class to be really immutable, we should not have a move constr
     or move assign operator. One way would be to declare them as deleted,
     but we do NOT want to do that.
     If we did that, the move cnstr/assign would still participate in OR,
     which would cause a compiler error in some cases, like when trying
     to move construct and object. So is there a better way? There is.
     We exploit the fact that this class has a user-declared copy constructor
     and copy assignment, so the compiler does not generate automatically
     a move constructor/move assignment, which means that only the copy
     constr/copy assign participate in overload resolution, which means we
     can achieve what we want by simply not declaring move cnstr/assign.

     See this for a full detailed explanation:
     https://blog.knatten.org/2021/10/15/the-difference-between-no-move-constructor-and-a-deleted-move-constructor/
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
