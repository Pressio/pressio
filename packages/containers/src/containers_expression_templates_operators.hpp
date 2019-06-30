
#ifndef CONTAINERS_EXPRESSION_TEMPLATES_OPERATORS_HPP_
#define CONTAINERS_EXPRESSION_TEMPLATES_OPERATORS_HPP_

namespace rompp{ namespace containers{ namespace exprtemplates{

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

}}} // end of namespace rompp::containers::exprtemplates
#endif
