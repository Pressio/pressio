
#ifndef ROMPP_MPL_PUBLICLY_INHERITS_FROM_HPP_
#define ROMPP_MPL_PUBLICLY_INHERITS_FROM_HPP_

namespace rompp{ namespace mpl{

template<typename T, typename base_t>
struct publicly_inherits_from : std::is_base_of<base_t,T>{};

}} // end namespace rompp::mpl
#endif
