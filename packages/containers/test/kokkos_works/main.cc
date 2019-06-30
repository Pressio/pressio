
#include <Kokkos_Core.hpp>
#include <cstdio>
#include <typeinfo>


struct hello_world {
  // If a functor has an "execution_space" (or "execution_space", for
  // backwards compatibility) public typedef, parallel_* will only run
  // the functor in that execution space.  That's a good way to mark a
  // functor as specific to an execution space.  If the functor lacks
  // this typedef, parallel_for will run it in the default execution
  // space, unless you tell it otherwise (that's an advanced topic;
  // see "execution policies").

  // The functor's operator() defines the loop body.  It takes an
  // integer argument which is the parallel for loop index.  Other
  // arguments are possible; see the "hierarchical parallelism" part
  // of the tutorial.
  //
  // The operator() method must be const, and must be marked with the
  // KOKKOS_INLINE_FUNCTION macro.  If building with CUDA, this macro
  // will mark your method as suitable for running on the CUDA device
  // (as well as on the host).  If not building with CUDA, the macro
  // is unnecessary but harmless.
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    printf ("Hello from i = %i\n", i);
  }
};

int main (int argc, char* argv[]) {
  // You must call initialize() before you may call Kokkos.
  //
  // With no arguments, this initializes the default execution space
  // (and potentially its host execution space) with default
  // parameters.  You may also pass in argc and argv, analogously to
  // MPI_Init().  It reads and removes command-line arguments that
  // start with "--kokkos-".
  Kokkos::initialize (argc, argv);

  // Print the name of Kokkos' default execution space.  We're using
  // typeid here, so the name might get a bit mangled by the linker,
  // but you should still be able to figure out what it is.
  printf ("Hello World on Kokkos execution space %s\n",
          typeid (Kokkos::DefaultExecutionSpace).name ());

  // Run the above functor on the default Kokkos execution space in
  // parallel, with a parallel for loop count of 15.
  //
  // The Kokkos::DefaultExecutionSpace typedef gives the default
  // execution space.  Depending on how Kokkos was configured, this
  // could be OpenMP, Threads, Cuda, Serial, or even some other
  // execution space.
  //
  // The following line of code would look like this in OpenMP:
  //
  // #pragma omp parallel for
  // for (int i = 0; i < 15; ++i) {
  //   printf ("Hello from i = %i\n", i);
  // }
  //
  // You may notice that the printed numbers do not print out in
  // order.  Parallel for loops may execute in any order.
  Kokkos::parallel_for ("HelloWorld",15, hello_world ());

  // You must call finalize() after you are done using Kokkos.
  Kokkos::finalize ();

  printf ("PASSED");
}

