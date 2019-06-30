
#ifndef ROMPPAPPS_BURGERS1D_EPETRA_MASKED_HPP_
#define ROMPPAPPS_BURGERS1D_EPETRA_MASKED_HPP_

#include "apps_burgers1d_epetra.hpp"

#ifdef HAVE_TRILINOS
#include <Epetra_Import.h>

namespace rompp{ namespace apps{

class Burgers1dEpetraMasked : public Burgers1dEpetra{
  using base_t = Burgers1dEpetra;
  using importer_t = Epetra_Import;

public:
  Burgers1dEpetraMasked(std::vector<scalar_type> params,
			int Ncell, Epetra_MpiComm * comm)
    : base_t(params, Ncell, comm){}

  ~Burgers1dEpetraMasked() = default;

public:
  void setup(){
    base_t::setup();
    // create a map to mimic the mask
    createMaskMap();
    //    maskMap_->Print(std::cout);
    importer_ = std::make_shared<importer_t>(*maskMap_, *dataMap_);
  };

  template <typename T>
  void applyMask(const T & src, T & dest, double t) const{
    dest.Import(src, *importer_, Insert);
  }

  residual_type applyMask(const residual_type & src, double t) const{
    residual_type dest(*maskMap_);
    dest.Import(src, *importer_, Insert);
    return dest;
  }

  Epetra_MultiVector applyMask(const Epetra_MultiVector & src,
			       double t) const{
    Epetra_MultiVector dest(*maskMap_, src.NumVectors());
    dest.Import(src, *importer_, Insert);
    return dest;
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

}} //namespace rompp::apps
#endif
#endif
