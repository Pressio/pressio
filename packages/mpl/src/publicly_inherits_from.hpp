
#ifndef PRESSIO_MPL_PUBLICLY_INHERITS_FROM_HPP_
#define PRESSIO_MPL_PUBLICLY_INHERITS_FROM_HPP_

namespace pressio{ namespace mpl{

template<typename T, typename base_t>
struct publicly_inherits_from : std::is_base_of<base_t,T>{};

}} // end namespace pressio::mpl
#endif
