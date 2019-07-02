
#ifndef QR_ALGORITHMS_TAGS_HPP_
#define QR_ALGORITHMS_TAGS_HPP_

namespace pressio{ namespace qr{

/* Note to devs: if you add a method here, and code it
 * remember to update the corresponding static_assert
 * message in the qr_traits.hpp */

struct ModifiedGramSchmidt{};
struct Householder{};

#if defined HAVE_TRILINOS
struct TSQR{};
#endif

}} // end namespace pressio::qr
#endif
