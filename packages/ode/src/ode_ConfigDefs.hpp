
#ifndef ODE_CONFIGDEFS_HPP_
#define ODE_CONFIGDEFS_HPP_

#include "core_ConfigDefs.hpp"
#include "meta/core_meta_basic.hpp"
//#include "core_crtp_helper.hpp"

namespace ode{
namespace details {

    template<typename T, typename enable = void>
    struct traits;
    // traits<const T> == traits<T> 
    // For example, traits<Vector<const T> > != traits<Vector<T> >, but
    //              traits<const Vector<T> > == traits<Vector<T> >
    template<typename T> 
    struct traits<const T> : traits<T> {};

    using time_type = double;

} // end namespace details
} // end ode

//namespace timeIntegrator{
// namespace defaultTypes {
// // namespace DefaultTypes {
// //! Default value of Scalar template parameter.
// using scalar_t = ::core::defaultTypes::scalar_t;
 // //! Default value of LocalOrdinal template parameter.
// //typedef int local_ordinal_t;
// /// \typedef global_ordinal_type
// //typedef int global_ordinal_t;
// } // namespace defaultTypes
//} // end of timeIntegrator namespace


#endif
