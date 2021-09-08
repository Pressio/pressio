
#ifndef TPQL_ROM_DATATYPES_HPP_
#define TPQL_ROM_DATATYPES_HPP_

namespace pressio{ namespace rom{ namespace experimental{

template<typename scalar_t>
struct RomDataTypeEigen
{
  using vector_t = Eigen::Matrix<scalar_t,-1,1>;
  using dense_matrix_t = Eigen::Matrix<scalar_t,-1,-1>;
  using dense_hessian_t = std::vector< std::vector< std::vector< scalar_t> > >; 
};


template<typename scalar_t>
struct RomDataTypeKokkos
{
  using execution_space = Kokkos::DefaultExecutionSpace;
  using kll   = Kokkos::LayoutLeft;
  using vector_t    = Kokkos::View<scalar_t*, kll, execution_space>;
  using dense_matrix_t    = Kokkos::View<scalar_t**, kll, execution_space>;
  using dense_hessian_t    = Kokkos::View<scalar_t***, kll, execution_space>;
};

}}}//end pressio::rom::experimental

#endif
