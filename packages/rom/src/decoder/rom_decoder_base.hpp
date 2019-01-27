
#ifndef ROM_DECODER_BASE_HPP_
#define ROM_DECODER_BASE_HPP_

namespace rompp{ namespace rom{

template <typename derived_type, typename jac_matrix_type>
struct DecoderBase
  : private core::details::CrtpBase<
  DecoderBase<derived_type, jac_matrix_type>>{

  using this_t = DecoderBase<derived_type, jac_matrix_type>;

  template <typename operand_t, typename result_t>
  void applyMapping(const operand_t & operandObj,
		    result_t & result) const{
    this->underlying().applyMappingImpl(operandObj, result);
  }

  const jac_matrix_type & getJacobianRef() const{
    return this->underlying().getJacobianRefImpl();
  }

private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  DecoderBase() = default;
  ~DecoderBase() = default;

};//end class

}}//end namespace rompp::rom
#endif
