
#ifndef CORE_CONFIGDEFS_HPP_
#define CORE_CONFIGDEFS_HPP_

#include "core_crtp_helper.hpp"
#include "core_config.h"
#include <type_traits>
#include "core_shared_traits.hpp"
#include "meta/core_meta_basic.hpp"
#include "core_default_types.hpp"
#include "core_static_constants.hpp"

#ifdef HAVE_TEUCHOS_TIMERS
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_StackedTimer.hpp"
#endif

namespace rompp{ namespace core{

namespace impl{

#ifdef DEBUG_PRINT
template <typename stream_t, typename scalar_t>
void setStreamPrecision(stream_t & ss){
  constexpr auto prec = std::is_same<scalar_t,
				     double>::value ? 15 : 6;
  ss << std::setprecision(prec);
}
#endif

struct empty{};

}//end namespace impl


namespace details {

template<typename T, typename enable = void>
struct traits : public
containers_shared_traits<void, void,
			 false, false, false,
			 WrappedPackageIdentifier::Undefined,
			 false, false>{};

template<typename T>
struct traits<const T> : traits<T> {};

} // end namespace details
//--------------------------------------------

namespace exprtemplates{

struct plus_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a+b) {
    return a + b;
  }
};

struct subtract_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a-b) {
    return a - b;
  }
};

struct times_{
  template <typename a_t, typename b_t>
  auto operator()(const a_t & a, const b_t & b) const
  -> decltype(a*b) {
    return a * b;
  }
};

} // end namespace exprtemplates
//--------------------------------------------


#ifdef HAVE_TEUCHOS_TIMERS
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
#endif

}} // end of namespace rompp::core
#endif
