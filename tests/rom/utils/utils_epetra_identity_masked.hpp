/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_epetra_identity_masked.hpp
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

#ifndef ROM_TESTS_BURGERS1D_EPETRA_IDENT_IDENT_MASKER_HPP_
#define ROM_TESTS_BURGERS1D_EPETRA_IDENT_IDENT_MASKER_HPP_

#include "pressio_apps.hpp"
#include <Epetra_Import.h>

namespace pressio{ namespace rom{ namespace test{ 

class Burgers1dEpetraIdentityMask //: public Burgers1dEpetra
{
public:
  using importer_t = Epetra_Import;
  // using base_t = Burgers1dEpetra;
  using dense_matrix_type = typename ::pressio::apps::Burgers1dEpetra::dense_matrix_type;
  using velocity_type = typename ::pressio::apps::Burgers1dEpetra::velocity_type;
  using state_type = typename ::pressio::apps::Burgers1dEpetra::state_type;
  using scalar_type = typename ::pressio::apps::Burgers1dEpetra::scalar_type;

  Burgers1dEpetraIdentityMask(Epetra_MpiComm * comm, 
    const Epetra_Map & dataMap)
  {
    createMaskMap(comm, dataMap);
    importer_ = std::make_shared<importer_t>(*maskMap_, dataMap);
    // this->setup();
  }

  ~Burgers1dEpetraIdentityMask() = default;

public:
  velocity_type createApplyMaskResult(const velocity_type & src) const
  {
    velocity_type dest(*maskMap_);
    dest.Import(src, *importer_, Insert);
    return dest;
  }

  dense_matrix_type createApplyMaskResult(const dense_matrix_type & src) const{
    dense_matrix_type dest(*maskMap_, src.NumVectors());
    dest.Import(src, *importer_, Insert);
    return dest;
  }

  template <typename T>
  void applyMask(const T & src, double time, T & dest) const{
    dest.Import(src, *importer_, Insert);
  }

private:
  void createMaskMap(Epetra_MpiComm * comm, const Epetra_Map & dataMap)
  {
    // get # of my elements for the full map
    auto myN0 = dataMap.NumMyElements();
    // get my global IDs
    auto myGID = dataMap.MyGlobalElements();

    // pick all elements
    std::vector<int> myGIDnc;
    for (decltype(myN0) i=0; i<myN0; i++) {
	   myGIDnc.emplace_back(myGID[i]);
    }
    maskMap_ = std::make_shared<Epetra_Map>(-1, myGIDnc.size(),
					    myGIDnc.data(), 0,
					    *comm);
    //maskMap_->Print(std::cout);
  };

private:
  template<typename T> using rcp = std::shared_ptr<T>;

  rcp<Epetra_Map> maskMap_;
  std::vector<int> myMaskGel_;
  rcp<importer_t> importer_;
};

}}}
#endif
