
#ifdef HAVE_TEUCHOS_TIMERS
#ifndef CORE_TEUCHOS_PERFORMANCE_MONITOR_HPP_
#define CORE_TEUCHOS_PERFORMANCE_MONITOR_HPP_

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_StackedTimer.hpp"

namespace rompp{ namespace core{

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
    timer->report(outp);
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

}} // end of namespace rompp::core
#endif
#endif
