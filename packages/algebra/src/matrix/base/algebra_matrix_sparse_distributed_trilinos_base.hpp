
#ifndef ALGEBRA_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_BASE_HPP_
#define ALGEBRA_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_BASE_HPP_

#include "../algebra_matrix_traits.hpp"

namespace rompp{
namespace algebra{

template<typename derived_type>
class MatrixSparseDistributedTrilinosBase
  : private utils::details::CrtpBase<
  MatrixSparseDistributedTrilinosBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==0,
  "OOPS: non-distributed matrix inheriting from sparse distributed trilinos!");

  using traits_t = details::traits<derived_type>;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using range_map_t = typename traits_t::range_map_t;
  using domain_map_t = typename traits_t::domain_map_t;

public:

  bool isFillingCompleted() const{
    return this->underlying().isFillingCompletedImpl();}

  void fillingIsCompleted(){
    this->underlying().fillingIsCompletedImpl();}

  void fillingIsCompleted(domain_map_t const & dmap,
			  range_map_t const & rmap){
    this->underlying().fillingIsCompletedImpl(dmap, rmap);}

  range_map_t const & getRangeDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getRangeDataMapImpl();}

  domain_map_t const & getDomainDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getDomainDataMapImpl();}

  row_map_t const & getRowDataMap() const{
    return this->underlying().getRowDataMapImpl();}

  col_map_t const & getColDataMap() const{
    return this->underlying().getColDataMapImpl();}

  bool hasSameRangeDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRangeDataMapAsImpl(other);}

  bool hasSameDomainDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameDomainDataMapAsImpl(other);}

  bool hasSameRowDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRowDataMapAsImpl(other);}

  bool hasSameColDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameColDataMapAsImpl(other);}

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<
    MatrixSparseDistributedTrilinosBase<derived_type>>;

  MatrixSparseDistributedTrilinosBase() = default;
  ~MatrixSparseDistributedTrilinosBase() = default;

};//end class

} // end namespace algebra
}//end namespace rompp
#endif
