
#ifndef ROM_DECODER_BASE_HPP_
#define ROM_DECODER_BASE_HPP_

#include "../rom_ConfigDefs.hpp"

namespace pressio{ namespace rom{

template <typename derived_type, typename jac_matrix_type>
struct DecoderBase
  // : private utils::details::CrtpBase<
  // DecoderBase<derived_type, jac_matrix_type>>
{
  using this_t = DecoderBase<derived_type, jac_matrix_type>;

  template <typename operand_t, typename result_t>
  void applyMapping(const operand_t & operandObj,
		    result_t & result) const  {
    static_cast<const derived_type &>(*this).applyMappingImpl(operandObj, result);
  }

  const jac_matrix_type & getReferenceToJacobian() const {
    return static_cast<const derived_type &>(*this).getReferenceToJacobianImpl();
  }

  DecoderBase() = default;
  ~DecoderBase() = default;

// private:
//   /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
//   template<typename DummyType> struct dummy{using type = DummyType;};
//   friend typename dummy<derived_type>::type;
//   friend utils::details::CrtpBase<this_t>;
};//end class

}}//end namespace pressio::rom
#endif
