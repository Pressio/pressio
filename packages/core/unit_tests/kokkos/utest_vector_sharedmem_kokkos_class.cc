
#include <gtest/gtest.h>
#include "CORE_VECTOR"

template <typename T>
struct InitView {
  T a;

  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  InitView (T a_) :
    a (a_){}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    // Acesss the View just like a Fortran array.  The layout depends
    // on the View's memory space, so don't rely on the View's
    // physical memory layout unless you know what you're doing.
    a[i] = 2.0*i;
  }
};

// template<typename T>
// struct ReduceFunctor {
//   T a;

//   // Constructor takes View by "value"; this does a shallow copy.
//   ReduceFunctor (T a_) : a (a_) {}

//   // If you write a functor to do a reduction, you must specify the
//   // type of the reduction result via a public 'value_type' typedef.
//   typedef double value_type;

//   KOKKOS_INLINE_FUNCTION
//   void operator() (int i, double &lsum) const {
//     lsum += a(i,0)*a(i,1)/(a(i,2)+0.1);
//   }
// };

// struct hello_world {
//   KOKKOS_INLINE_FUNCTION
//   void operator() (const int i) const {
//     printf ("Hello from i = %i\n", i);
//   }
// };


TEST(core_vector_sharedmem_kokkos_class, Constructor)
{
  // kokkos initialize and finalize already set from environment, see CMakeList
  const int N = 10;

  using view_type = Kokkos::View<double*>;
  view_type a ("A", N);
  //Kokkos::parallel_for ("HelloWorld",15, hello_world());
  
  using myvec_t = core::Vector<view_type>;
  static_assert( core::details::traits<myvec_t>::isStatic == 0, "" );
  myvec_t g(a);
  
  Kokkos::parallel_for (N, InitView<view_type>( g.dataCp() ));
  // double sum = 0;
  // Kokkos::parallel_reduce (N, ReduceFunctor<view_type>(a), sum);
  // printf ("Result: %f\n", sum);


  using view_type2 = Kokkos::View<double[11]>;
  using myvec_t2 = core::Vector<view_type2>;
  static_assert( core::details::traits<myvec_t2>::isStatic == 1, "" );
  
}
