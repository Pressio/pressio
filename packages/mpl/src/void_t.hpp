
#ifndef ROMPP_MPL_VOID_T_HPP_
#define ROMPP_MPL_VOID_T_HPP_

namespace rompp { namespace mpl{

// A void_t implementation that works with gcc-4.9 (workaround for bug 64395
// From: http://stackoverflow.com/questions/35753920/why-does-the-void-t-detection-idiom-not-work-with-gcc-4-9
namespace _void_t_impl {

template <class... >
struct make_void { using type = void; };

} // end namepace _void_t_impl

template <class... T>
using void_t = typename _void_t_impl::make_void<T...>::type;

}} // end namespace 
#endif 
