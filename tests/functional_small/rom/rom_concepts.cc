
#include <gtest/gtest.h>
#include "pressio/rom_decoder.hpp"

namespace
{
  struct ValidFomState{
    ValidFomState(const ValidFomState &) = default;
  };

  struct InvalidFomState{
    InvalidFomState(const InvalidFomState &) = delete;
  };


template <class MatrixType, class FomStateType>
struct ValidDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = MatrixType;
  using fom_state_type = FomStateType;

public:
  template <typename RomStateType>
  void applyMapping(const RomStateType & romState, 
                    fom_state_type & result) const;

  const jacobian_type & jacobianCRef() const; 

  template <typename RomStateType>
  void updateJacobian(const RomStateType &);
};//end

}//end namespace anonim

TEST(rom, concepts_fom_state)
{
  namespace prom = pressio::rom;
  static_assert(prom::fom_state<ValidFomState>::value, "");
  static_assert(!prom::fom_state<InvalidFomState>::value, "");
}

TEST(rom, concepts_decoder)
{
  namespace prom = pressio::rom;

  struct AMatrixClass{};
  struct AVectorClass{};
  using decoder_t = ValidDecoder<AMatrixClass, AVectorClass>;
  using rom_state_t = std::vector<double>;
  static_assert(prom::decoder<decoder_t, rom_state_t>::value, "");
}