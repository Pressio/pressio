/*
//@HEADER
// ************************************************************************
//
// utils_teuchos_performance_monitor.hpp
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

#ifndef UTILS_UTILS_TEUCHOS_PERFORMANCE_MONITOR_HPP_
#define UTILS_UTILS_TEUCHOS_PERFORMANCE_MONITOR_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_StackedTimer.hpp"

namespace pressio{ namespace utils{

struct TeuchosPerformanceMonitor{

  using base_comm_t = Teuchos::Comm<int>;
  using rcp_comm_t  = Teuchos::RCP<const base_comm_t>;
  using ser_comm_t  = Teuchos::SerialComm<int>;
  using mpi_comm_t  = Teuchos::MpiComm<int>;
  using def_comm_t  = Teuchos::DefaultComm<int>;

  static
  void stackedTimersReportSerial(std::ostream & outp = std::cout){

    Teuchos::StackedTimer::OutputOptions options;
    options.output_histogram=true;
    options.output_fraction=true;
    options.output_minmax=true;

    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->report(outp, Teuchos::rcp(new ser_comm_t()), options);
  }

  static
  void stackedTimersReportMPI(std::ostream & outp = std::cout,
			      rcp_comm_t comm = def_comm_t::getComm()){
    /*Teuchos::(new mpi_comm_t(MPI_COMM_WORLD))){*/

    Teuchos::StackedTimer::OutputOptions options;
    //options.output_histogram=true;
    options.output_fraction=true;
    options.output_minmax=true;

    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->report(outp, comm, options);

    // auto mystream = Teuchos::rcp(new std::ofstream("timing.log",
    // 						 std::ios::out));
    // Teuchos::FancyOStream timing_stream(mystream);
    // timer->report(timing_stream, comm, options);
  }
};

}} // end of namespace pressio::utils
#endif  // UTILS_UTILS_TEUCHOS_PERFORMANCE_MONITOR_HPP_
