
#ifndef PRESSIOAPPS_BURGERS1D_EPETRA_REDUCED_NO_MASK_HPP_
#define PRESSIOAPPS_BURGERS1D_EPETRA_REDUCED_NO_MASK_HPP_

#include "apps_burgers1d_epetra.hpp"

#ifdef HAVE_TRILINOS
#include <Epetra_Import.h>

namespace pressio{ namespace apps{

class Burgers1dEpetraReducedNoMask : public Burgers1dEpetra{
  using base_t	   = Burgers1dEpetra;
  using importer_t = Epetra_Import;
  using MV_t	   = Epetra_MultiVector;

public:
  Burgers1dEpetraReducedNoMask(std::vector<scalar_type> params,
			int Ncell, Epetra_MpiComm * comm)
    : base_t(params, Ncell, comm){}

  ~Burgers1dEpetraReducedNoMask() = default;

public:
  void setup(){
    base_t::setup();
    // create a map to mimic the mask
    createMaskMap();
    //    maskMap_->Print(std::cout);
    importer_ = std::make_shared<importer_t>(*maskMap_, *dataMap_);
  };

  void velocity(const state_type & u,
		const scalar_type t, 
    velocity_type & rhs) const
  {
    velocity_type R(*dataMap_);
    base_t::velocity(u, t, R);
    rhs.Import(R, *importer_, Insert);
  }

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const{
    velocity_type R(*dataMap_);
    base_t::velocity(u, t, R);
    velocity_type dest(*maskMap_);
    dest.Import(R, *importer_, Insert);
    return dest;
  }

  // A = Jac * B
  void applyJacobian(const state_type & y,
		     const MV_t & B,
		     scalar_type t,
         MV_t & A) const{
    MV_t Cfull( Jac_->RangeMap(), B.NumVectors() );
    base_t::applyJacobian(y, B, t, Cfull);
    A.Import(Cfull, *importer_, Insert);
  }

  // return Jac * B
  MV_t applyJacobian(const state_type & y,
		     const MV_t & B,
		     scalar_type t) const{
    MV_t Cfull( Jac_->RangeMap(), B.NumVectors() );
    base_t::applyJacobian(y, B, t, Cfull);
    MV_t C(*maskMap_, B.NumVectors());
    C.Import(Cfull, *importer_, Insert);
    return C;
  }

private:
  void createMaskMap(){
    // get # of my elements for the full map
    auto myN0 = dataMap_->NumMyElements();
    // get my global IDs
    auto myGID = dataMap_->MyGlobalElements();

    // pick some elements
    std::vector<int> myGIDnc;
    for (decltype(myN0) i=0; i<myN0; i++) {
      if ( i<=8 or (i>=11 and i<=18) or i>=21 )
	myGIDnc.emplace_back(myGID[i]);
    }

    maskMap_ = std::make_shared<Epetra_Map>(-1, myGIDnc.size(),
					    myGIDnc.data(), 0,
					    *comm_);
    maskMap_->Print(std::cout);
  };

private:
  rcp<Epetra_Map> maskMap_;
  std::vector<int> myMaskGel_;
  rcp<importer_t> importer_;
};

}} //namespace pressio::apps
#endif
#endif
