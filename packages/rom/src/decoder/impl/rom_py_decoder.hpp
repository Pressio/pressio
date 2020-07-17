
#ifndef ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_
#define ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  typename matrix_type,
  typename rom_state_type,
  typename fom_state_type
  >
struct PyDecoder
{
  using jacobian_type  = matrix_type;
  enum mappingKind{ Null, Linear, Custom };

private:
  using scalar_t	  = typename ::pressio::containers::details::traits<rom_state_type>::scalar_t;
  using jacobian_native_t = typename ::pressio::containers::details::traits<jacobian_type>::wrapped_t;
  using fom_native_t	  = typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t;

  matrix_type mappingJacobian_ = {};
  mappingKind kind_ = {};
  pybind11::object customMapper_ = {};

public:
  PyDecoder() = delete;

  PyDecoder(const jacobian_native_t & matIn)
    : mappingJacobian_(matIn), kind_(mappingKind::Linear)
  {}

  PyDecoder(pybind11::object customMapper)
    //note here that we view the native object, we don't deep copy it
    // so if the mapping jacobian changes on the python side, it is also reflected here
    : mappingJacobian_(customMapper.attr("getJacobian")(), ::pressio::view()),
      kind_(mappingKind::Custom),
      customMapper_(customMapper)
  {}

  // applyMapping is templated because operand_t can be rom_state_type but
  // can also be an expression based on rom_state_type (e.g. for WLS)
  template <typename operand_t, typename fom_state_t = fom_state_type>
  void applyMapping(const operand_t & operand, fom_state_t & result) const
  {
    if (kind_ == mappingKind::Linear){
      constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
      constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
      ::pressio::ops::product(::pressio::nontranspose(), one, mappingJacobian_, operand, zero, result);
    }
    else if(kind_ == mappingKind::Custom)
    {
      customMapper_.attr("applyMapping")(*operand.data(), *result.data());
    }
    else
      throw std::runtime_error("Invalid mapping kind enum");
  }

  const jacobian_type & getReferenceToJacobian() const{
    return mappingJacobian_;
  }
};//end

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_
