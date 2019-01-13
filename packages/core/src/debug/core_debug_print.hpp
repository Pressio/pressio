
#ifndef CORE_DEBUG_PRINT_HPP_
#define CORE_DEBUG_PRINT_HPP_

#include "../core_ConfigDefs.hpp"
#include <iostream>

namespace rompp{ namespace core{ namespace debug{

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

}//end namepsace impl


template <typename ... Args>
void print(Args... args){
  int myRank = 0;
#if defined HAVE_MPI
  int flag = 0; MPI_Initialized( &flag );
  if (flag==1) MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  if(myRank==0)
    impl::print(std::cout, std::forward<Args>(args)... );
}

}}}//end namespace rompp::core::debug
#endif
