
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
void print(Args &&... args){
  if(impl::myRank()==0)
    impl::print(std::cout, std::forward<Args>(args)... );
}

template <typename T, typename ... Args>
void printCoreWrapper(const T & object, Args && ...args)
{
  using traits = core::details::traits<T>;

  if (traits::is_shared_mem == 1){
    if(impl::myRank()==0)
      object.print(std::forward<Args>(args)...);
  }
  else
    object.print(std::forward<Args>(args)...);
}


}}}//end namespace rompp::core::debug
#endif
