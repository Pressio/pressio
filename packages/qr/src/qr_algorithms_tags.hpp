
#ifndef QR_ALGORITHMS_TAGS_HPP_
#define QR_ALGORITHMS_TAGS_HPP_

namespace rompp{ namespace qr{

struct ModifiedGramSchmidt{};
struct Householder{};

#if defined HAVE_TRILINOS
struct TSQR{};
#endif

}} // end namespace rompp::qr
#endif
