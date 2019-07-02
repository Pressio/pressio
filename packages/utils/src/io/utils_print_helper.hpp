
#ifndef UTILS_PRINT_HELPER_HPP_
#define UTILS_PRINT_HELPER_HPP_

#include "../utils_ConfigDefs.hpp"
#include "utils_colorize_print.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace pressio{ namespace utils{ namespace io{

namespace impl{

template <typename head_t>
void print(std::ostream & ss, head_t && head){
  ss << std::forward<head_t>(head);
}

template <typename head_t, typename ... tail_t>
void print(std::ostream & ss, head_t && head, tail_t && ... args){
  ss << std::forward<head_t>(head) << " ";
  print(ss, std::forward<tail_t>(args)...);
}

template <typename T = int>
T myRank(){
  T rank = 0;
#if defined HAVE_MPI
  int flag = 0; MPI_Initialized( &flag );
  if (flag==1) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return static_cast<T>(rank);
}

}//end namepsace impl


template <typename ... Args>
void print_stdout(Args &&... args){
  if(impl::myRank()==0)
    impl::print(std::cout, std::forward<Args>(args)... );
}

// conditional on DEBUG_PRINT being on
template <typename ... Args>
void print_stdout_if_dp(Args &&... args){
#if defined DEBUG_PRINT
  print_stdout(std::forward<Args>(args)...);
#endif
}


// template <typename T, typename ... Args>
// void print_containers_wrapper(const T & object, Args && ...args)
// {
//   using traits = containers::details::traits<T>;

//   if (traits::is_shared_mem == true){
//     if(impl::myRank()==0)
//       object.print(std::forward<Args>(args)...);
//   }
//   else
//     object.print(std::forward<Args>(args)...);
// }


}}}//end namespace pressio::utils::io
#endif
