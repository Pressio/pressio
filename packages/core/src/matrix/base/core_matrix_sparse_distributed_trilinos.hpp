
#ifndef CORE_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_HPP_
#define CORE_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class MatrixSparseDistributedTrilinos
  : private core::details::CrtpBase<
  MatrixSparseDistributedTrilinos<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==1,
		 "OOPS: non-distributed matrix inheriting \
from sparse distributed trilinos!");

private:
  using traits_t = details::traits<derived_type>;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using range_map_t = typename traits_t::range_map_t;
  using domain_map_t = typename traits_t::domain_map_t;
  // using crs_graph_t = typename traits_t::crs_graph_t;

public:
  
  bool isFillingCompleted() const{
    return this->underlying().isFillingCompletedImpl();
  }

  void fillingIsCompleted(){
    this->underlying().fillingIsCompletedImpl();
  }

  range_map_t const & getRangeDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getRangeDataMapImpl();
  }
 
  domain_map_t const & getDomainDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getDomainDataMapImpl();
  }

  row_map_t const & getRowDataMap() const{
    return this->underlying().getRowDataMapImpl();
  }
 
  col_map_t const & getColDataMap() const{
    return this->underlying().getColDataMapImpl();
  }

  bool hasSameRangeDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRangeDataMapAsImpl(other);
  }

  bool hasSameDomainDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameDomainDataMapAsImpl(other);
  }

  bool hasSameRowDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRowDataMapAsImpl(other);
  }

  bool hasSameColDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameColDataMapAsImpl(other);
  }
  
  // crs_graph_t const & getCrsGraph() const{
  //   return this->underlying().getCrsGraphImpl();
  // }

private:
  friend derived_type;
  friend core::details::CrtpBase<
    MatrixSparseDistributedTrilinos<derived_type>>;

  MatrixSparseDistributedTrilinos() = default;
  ~MatrixSparseDistributedTrilinos() = default;
  
};//end class
  
} // end namespace core
#endif
