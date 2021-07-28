/*
//@HEADER
// ************************************************************************
//
// utils_print_helper.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef UTILS_IO_UTILS_PRINT_HELPER_HPP_
#define UTILS_IO_UTILS_PRINT_HELPER_HPP_

#ifdef PRESSIO_ENABLE_TPL_MPI
#include <mpi.h>
#endif
#include <iostream>

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
#if defined PRESSIO_ENABLE_TPL_MPI
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

// // conditional on PRESSIO_ENABLE_DEBUG_PRINT being on
// template <typename ... Args>
// void print_stdout_if_dp(Args &&... args){
// #if defined PRESSIO_ENABLE_DEBUG_PRINT
//   print_stdout(std::forward<Args>(args)...);
// #endif
// }

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
#endif  // UTILS_IO_UTILS_PRINT_HELPER_HPP_
