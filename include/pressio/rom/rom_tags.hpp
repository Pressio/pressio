
#ifndef ROM_TAGS_HPP_
#define ROM_TAGS_HPP_

namespace pressio{ namespace rom{

namespace impl{
struct FomStepTag{};
struct RomStepTag{};
struct TransitionToRomAndDoStepTag{};
struct TransitionToFomAndDoStepTag{};
}

constexpr impl::FomStepTag FomStep{};
constexpr impl::RomStepTag RomStep{};
constexpr impl::TransitionToRomAndDoStepTag TransitionToRomAndDoStep{};
constexpr impl::TransitionToFomAndDoStepTag TransitionToFomAndDoStep{};

}} // end pressio::rom
#endif
