/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_discrete_time_default_system.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ROM_IMPL_LSPG_UNSTEADY_RECONSTRUCTOR_HPP_
#define ROM_IMPL_LSPG_UNSTEADY_RECONSTRUCTOR_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<typename reduced_state_type>
auto read(std::string const & filein, std::size_t ext)
{
  using scalar_type = scalar_trait_t<reduced_state_type>;

  std::ifstream reduced_state_file(filein);
  assert(reduced_state_file.good());

  std::string line;
  std::vector<scalar_type> times;
  std::vector<reduced_state_type> reduced_states;
  while (std::getline(reduced_state_file, line)){
    std::istringstream stm(line);
    scalar_type time;
    reduced_state_type reduced_state(ext);
    scalar_type reduced_state_i;
    stm >> time;
    times.push_back(time);

    int c=0;
    while (stm >> reduced_state_i){
      reduced_state[c++] = reduced_state_i;
    }
    reduced_states.push_back(reduced_state);
  }
  reduced_state_file.close();
  return std::make_tuple(times, reduced_states);
}

template<typename ViewT>
auto _rank1_view_to_stdvector(ViewT view)
{
  using sc_t = std::remove_cv_t<scalar_trait_t<ViewT>>;
  std::vector<sc_t> v(view.extent(0));
  for (std::size_t k=0; k<view.extent(0); ++k){
    v[k] = view(k,0);
  }
  return v;
}

template<typename ViewT>
auto _rank2_view_to_stdvector(ViewT view)
{
  using sc_t = std::remove_cv_t<scalar_trait_t<ViewT>>;
  std::vector<sc_t> v(view.extent(0)*view.extent(1));
  std::size_t k=0;
  for (std::size_t j=0; j<view.extent(1); ++j){
    for (std::size_t i=0; i<view.extent(0); ++i){
      v[k++] = view(i,j);
    }
  }
  return v;
}

template<typename sc_t, typename size_t>
void write_vector_to_binary(const std::string filename,
			    const sc_t * A,
			    size_t n)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  out.write((char*) A, n*sizeof(sc_t) );
  out.close();
}

template <
  std::size_t n,
  class TrialSubspaceType,
  class FomSystemType
  >
class LspgReconstructor{
  using state_type = typename TrialSubspaceType::reduced_state_type;
  using scalar_type = scalar_trait_t<state_type>;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;

public:
  LspgReconstructor(const TrialSubspaceType & trialSubspace,
		    const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace), fomSystem_(fomSystem){}

  void execute(std::string const & filename) const
  {
    static_assert(n==2);

    const auto & trialSub = trialSubspace_.get();
    auto y_np1 = trialSub.createFullState();
    auto y_n   = trialSub.createFullState();
    const auto & phi = trialSub.basisOfTranslatedSpace();

    auto R = fomSystem_.get().createDiscreteTimeResidual();
    auto JPhi = fomSystem_.get().createResultOfDiscreteTimeJacobianActionOn(phi);
    assert( *R.getMap() == *JPhi.getMap() );

    const auto myRank = R.getMap()->getComm()->getRank();
    std::ofstream outMap("row_map.txt");
    auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(outMap));//std::cout));
    R.getMap()->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
    outMap.close();

    const std::size_t numModes = trialSubspace_.get().dimension();
    auto [times, redStates] = read<state_type>(filename, numModes);

    trialSub.mapFromReducedState(redStates[0], y_n);
    for (std::size_t i = 1; i < times.size(); i++){
      trialSub.mapFromReducedState(redStates[i], y_np1);
      fomSystem_.get().discreteTimeResidualAndJacobianAction(i,
							     times[i],
							     times[i] - times[i - 1],
							     R, phi, &JPhi,
							     y_np1, y_n);
      auto r_view = R.getLocalViewHost(Tpetra::Access::ReadOnly);
      auto r_stdv = _rank1_view_to_stdvector(r_view);
      const std::string R_f = "residual_rank_" + std::to_string(myRank) + "_" + std::to_string(i) + ".bin";
      write_vector_to_binary(R_f, r_stdv.data(), r_stdv.size());

      auto jphi_view = JPhi.getLocalViewHost(Tpetra::Access::ReadOnly);
      auto jphi_stdv = _rank2_view_to_stdvector(jphi_view);
      std::string Jphi_f = "jacobian_action_rank_" + std::to_string(myRank) + "_" + std::to_string(i) + ".bin";
      write_vector_to_binary(Jphi_f, jphi_stdv.data(), jphi_stdv.size());

      pressio::ops::deep_copy(y_n, y_np1);
    }
  }
};

}}}
#endif
