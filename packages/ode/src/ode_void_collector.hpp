
#ifndef ODE_VOID_COLLECTOR_HPP_
#define ODE_VOID_COLLECTOR_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{

struct voidCollector
{
    template< class state_t , class time_t >
    void operator()( size_t /*step*/,
		     time_t /*t*/,
		     const state_t & /* x */ )
    {}
};

} // end ode 
#endif

