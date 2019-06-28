
#include <gtest/gtest.h>
#include "ALGEBRA_VECTOR"

template <typename T>
struct InitView {
  T a;

  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  InitView (T a_) :
    a (a_)
  {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    // Acesss the View just like a Fortran array.  The layout depends
    // on the View's memory space, so don't rely on the View's
    // physical memory layout unless you know what you're doing.
    a(i,0) = 1.0*i;
    a(i,1) = 1.0*i*i;
    a(i,2) = 1.0*i*i*i;
  }
};

template<typename T>
struct ReduceFunctor {
  T a;

  // Constructor takes View by "value"; this does a shallow copy.
  ReduceFunctor (T a_) : a (a_) {}

  // If you write a functor to do a reduction, you must specify the
  // type of the reduction result via a public 'value_type' typedef.
  typedef double value_type;

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, double &lsum) const {
    lsum += a(i,0)*a(i,1)/(a(i,2)+0.1);
  }
};



TEST(algebra_vector_sharedmem_kokkos_class, Constructor)
{
  Kokkos::initialize();//argc, argv);
  const int N = 10;

  using view_type = Kokkos::View<double*[3]>;
  view_type a ("A", N);
  
  Kokkos::parallel_for (N, InitView<view_type>(a));
  double sum = 0;
  Kokkos::parallel_reduce (N, ReduceFunctor<view_type>(a), sum);
  printf ("Result: %f\n", sum);
  Kokkos::finalize ();
  
  Kokkos::finalize ();
  

  // // using eigvec_t1 = Eigen::Matrix<double, 4, 1>;
  // // using myvec_t1 = algebra::Vector<eigvec_t1>;
  // // myvec_t1 m333(5);
  // using eigvec_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  // using myvec_t = algebra::Vector<eigvec_t>;
  // using vecTrait = algebra::details::traits<myvec_t>;
  // ASSERT_TRUE(vecTrait::isEigen == 1);
 
  // myvec_t m_v2(5);
  // ASSERT_FALSE( m_v2.empty() );
  // ASSERT_TRUE( m_v2.size() == 5 );

  // // create an eigen-type vector
  // Eigen::Vector4d e_v1;
  // myvec_t m_v3(e_v1);
  // ASSERT_FALSE( m_v3.empty() );
  // ASSERT_TRUE( m_v3.size() == 4 );
}
