
#ifndef PRESSIO_ROM_LINEAR_SUBSPACE_HPP_
#define PRESSIO_ROM_LINEAR_SUBSPACE_HPP_

namespace pressio{ namespace rom{

template <class BasisType, class OffsetType>
class LinearAffineSubspace
{

public:
  using basis_type  = std::remove_cv_t<BasisType>;
  using offset_type = std::remove_cv_t<OffsetType>;

private:
  static_assert( !::pressio::mpl::all_of_t<std::is_pointer, basis_type, offset_type>::value,
		"template parameters cannot be pointers");
  static_assert( !::pressio::mpl::all_of_t<mpl::is_std_shared_ptr, basis_type, offset_type>::value,
		"std::unique_ptr are not valid template parameters");
  static_assert(std::is_copy_constructible< basis_type >::value,
		"BasisType must be copy constructible");
  static_assert(std::is_copy_constructible< offset_type >::value,
		"OffsetType must be copy constructible");

  // mandates
  static_assert(std::is_same<
		typename pressio::Traits< basis_type>::scalar_type,
		typename pressio::Traits< offset_type >::scalar_type >::value,
		"Mismatching scalar_type");

private:
  const basis_type basis_;
  const offset_type offset_;
  bool isAffine_ = {};

public:
  LinearAffineSubspace() = delete;

  LinearAffineSubspace(const basis_type & basis,
		       const offset_type & offset,
		       bool isAffine)
    : basis_(::pressio::ops::clone(basis)),
      offset_(::pressio::ops::clone(offset)),
      isAffine_(isAffine){}

  LinearAffineSubspace(basis_type && basis,
		       offset_type && offset,
		       bool isAffine)
    : basis_(std::move(basis)),
      offset_(std::move(offset)),
      isAffine_(isAffine){}

  LinearAffineSubspace(const basis_type & basis,
		       offset_type && offset,
		       bool isAffine)
    : basis_(::pressio::ops::clone(basis)),
      offset_(std::move(offset)),
      isAffine_(isAffine){}

  LinearAffineSubspace(basis_type && basis,
		       const offset_type & offset,
		       bool isAffine)
    : basis_(std::move(basis)),
      offset_(::pressio::ops::clone(offset)),
      isAffine_(isAffine){}

  const bool isAffine() const{ return isAffine_; }
  const basis_type & viewBasis() const{ return basis_; }
  const offset_type & viewOffset() const{ return offset_; }
};

}}
#endif
