
#ifndef CORE_EXPRESSION_TEMPLATES_OPERATORS_HPP_
#define CORE_EXPRESSION_TEMPLATES_OPERATORS_HPP_

namespace rompp{ namespace core{ namespace exprtemplates{

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

}}} // end of namespace rompp::core::exprtemplates
#endif
