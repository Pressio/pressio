
#ifndef ROMPPAPPS_BURGERS1D_TPETRA_HPP_
#define ROMPPAPPS_BURGERS1D_TPETRA_HPP_

#include "../../../CORE_ALL"
#include "../apps_ConfigDefs.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_FancyOStream.hpp>

namespace rompp{ namespace apps{

class Burgers1dTpetra{
protected:
  using map_t		= Tpetra::Map<>;
  using nativeVec	= Tpetra::Vector<>;
  using nativeMVec	= Tpetra::MultiVector<>;
  using jacobian_type	= Tpetra::CrsMatrix<>;
  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;

  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  template<typename T> using stdrcp = std::shared_ptr<T>;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= nativeVec;
  using residual_type	= state_type;

public:
  Burgers1dTpetra(std::vector<scalar_type> params,
		  int Ncell,
		  rcpcomm_t comm)
    : mu_{params}, Ncell_{Ncell}, comm_{comm}{}

  Burgers1dTpetra() = delete;
  ~Burgers1dTpetra() = default;

public:
  rcpmap_t getDataMap(){
    return dataMap_;
  };

  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));

    myRank_ =  comm_->getRank();
    totRanks_ =  comm_->getSize();

    // distribute cells
    dataMap_ = Teuchos::rcp(new map_t(Ncell_, 0, comm_));

    NumMyElem_ = dataMap_->getNodeNumElements();
    auto minGId = dataMap_->getMinGlobalIndex();
    myGel_.resize(NumMyElem_);
    std::iota(myGel_.begin(), myGel_.end(), minGId);
    dataMap_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    if (myRank_==1){
      for (auto it : myGel_)
	std::cout << it << std::endl;
    }

    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;

    // grid
    xGrid_ = std::make_shared<nativeVec>(dataMap_);
    auto xGridv = xGrid_->getDataNonConst();
    lo_t i=0;
    for (auto const & it : myGel_){
      xGridv[i] = dx_*it + dx_*0.5;
      i++;
    };
    xGrid_->describe(*wrappedCout_, Teuchos::VERB_EXTREME);

    // init condition
    U0_ = std::make_shared<nativeVec>(dataMap_);
    U0_->putScalar(1.0);

    U_ = std::make_shared<nativeVec>(dataMap_);
    U_->putScalar(1.0);

    Jac_ = std::make_shared<jacobian_type>(dataMap_, nonZrPerRow_);
    this->computeJacobianInsert(*U_, *Jac_, 0.0);
  };

  state_type const & getInitialState() const{
    return *U0_;
  };

  void residual(const state_type & u,
		residual_type & rhs,
		const scalar_type /* t */) const;

  residual_type residual(const state_type & u,
			 const scalar_type t) const{
    residual_type R(dataMap_);
    residual(u,R,t);
    return R;
  }

  // computes: A = Jac B where B is a multivector
  void applyJacobian(const state_type & y,
		     const nativeMVec & B,
		     nativeMVec & A,
		     scalar_type t) const{
    // assert( Jac_->NumGlobalCols() == B.GlobalLength() );
    // assert( A.GlobalLength() == Jac_->NumGlobalRows() );
    // assert( A.NumVectors() == B.NumVectors() );
    computeJacobianReplace(y, *Jac_, t);
    Jac_->apply(B, A);
    //A.Print(std::cout);
  }

  // computes: A = Jac B where B is a multivector
  nativeMVec applyJacobian(const state_type & y,
			   const nativeMVec & B,
			   scalar_type t) const{
    nativeMVec C( dataMap_, B.getNumVectors() );
    applyJacobian(y, B, C, t);
    return C;
  }

protected:
  void computeJacobianReplace(const state_type & u,
			      jacobian_type & jac,
			      const scalar_type /*t*/) const;

  void computeJacobianInsert(const state_type & u,
			    jacobian_type & jac,
			    const scalar_type /*t*/) const;

  scalar_type exchangeFlux(const state_type & u) const;

protected:
  std::vector<scalar_type> mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain
  const scalar_type xR_ = 100.0; // right side of domain
  int Ncell_{}; // # of cells
  scalar_type dx_{}; // cell size
  scalar_type dxInv_{}; // inv of cell size
  const int nonZrPerRow_ = 2;
  stdrcp<nativeVec> xGrid_{}; // mesh points coordinates

  Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
  rcpcomm_t comm_{};
  rcpmap_t dataMap_{};

  int myRank_{};
  int totRanks_{};
  lo_t NumMyElem_{};
  std::vector<go_t> myGel_{};

  mutable stdrcp<nativeVec> U_{}; // state vector
  mutable stdrcp<nativeVec> U0_{}; // initial state vector
  stdrcp<jacobian_type> Jac_{};

};//end class

}} //namespace rompp::apps
#endif
