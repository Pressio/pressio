
#ifdef HAVE_PYBIND11
#ifndef ROM_LINEAR_PY_DECODER_HPP_
#define ROM_LINEAR_PY_DECODER_HPP_

#include "rom_decoder_base.hpp"
#include "../rom_fwd.hpp"

namespace rompp{ namespace rom{

template <
  typename matrix_type,
  typename ops_t
  >
struct PyLinearDecoder<
  matrix_type, ops_t,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_cstyle_array_pybind11<matrix_type>::value
    >
  >
  : public DecoderBase<
  PyLinearDecoder<matrix_type, ops_t>, matrix_type
  >
{

  using this_t	    = PyLinearDecoder<matrix_type, ops_t>;
  using base_t	    = DecoderBase<this_t, matrix_type>;
  using jacobian_t  = matrix_type;
  using scalar_t    = double;
  static_assert( mpl::is_same<scalar_t,
		 typename matrix_type::value_type>::value,
		 "PyLinearDecoder: Scalar types don't match");

private:
  friend base_t;
  matrix_type phi_ = {};
  ops_t ops_ = {};

public:
  PyLinearDecoder() = delete;

  // since matrix is a pybind11::array_t, the following
  // uses view semantics, so NOT a deep copy. Hence, the object
  // is owened on the python side. This is fine, beucase this is the
  // only object who owns the basis vectors. if we want to do a deep copy,
  // then we have to init phi_ using the buffer of matIn like we do
  // in ode_storage for example
  PyLinearDecoder(const jacobian_t & matIn,
		  const ops_t ops)
    : phi_(matIn), ops_{ops}
  {
    std::cout << std::endl;
    std::cout << "Inside PyLinearDecoder " << std::endl;
    std::cout << "phi_ " << phi_.data() << std::endl;
    std::cout << std::endl;
  }

  ~PyLinearDecoder() = default;

protected:
  template <typename operand_t, typename result_t>
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const{
    ops_.attr("multiply2")(phi_, operandObj, resultObj);
  }

  const jacobian_t & getReferenceToJacobianImpl() const{
    return phi_;
  }

};//end class

}}//end namespace rompp::rom
#endif
#endif
