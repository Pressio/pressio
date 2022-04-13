
#include <gtest/gtest.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#include <KokkosBlas.hpp>
#include "pressio/utils.hpp"

template<class VecType>
void write_tpetra_vector_to_ascii_file(std::string fileName,
				       const VecType & tpVec)
{
  auto view = tpVec.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<view.extent(0); i++){
    file << std::setprecision(15) << view(i,0) << " \n";
  }
  file.close();
}

template<class DataType, class ...Prop>
void write_kokkos_view_to_ascii_file(std::string fileName,
				     const Kokkos::View<DataType, Prop...> kv)
{
  std::ofstream file; file.open(fileName);
  for (size_t i=0; i<kv.extent(0); i++){
    file << std::setprecision(15) << kv(i) << " \n";
  }
  file.close();
}

template<class T>
void read_integers_from_ascii_into_view(const std::string & fileName,
					Kokkos::View<T*> & dest)
{
  std::vector<int> v;
  std::ifstream source;
  source.open(fileName, std::ios_base::in);
  std::string line, colv;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    in >> colv;
    v.push_back(std::atoi(colv.c_str()));
  }
  source.close();

  Kokkos::resize(dest, v.size());
  for (int i=0; i<v.size(); ++i){
    dest(i) = v[i];
  }
}

template<class ViewType>
void read_phi_from_ascii_into_view(std::string fileName,
				   ViewType & M,
				   std::vector<int> targetRows)
{
  std::vector<std::vector<double>> A0;
  pressio::utils::read_ascii_matrix_stdvecvec(fileName, A0, M.extent(1));

  for (int i=0; i<targetRows.size(); ++i){
    for (int j=0; j<M.extent(1); ++j){
      M(i,j) = A0[targetRows[i]][j];
    }
  }
}

auto create_stdvector_and_fill_from_ascii(const std::string & fileName)
{
  std::vector<int> v;
  std::ifstream source;
  source.open(fileName, std::ios_base::in);
  std::string line, colv;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    in >> colv;
    v.push_back(std::atoi(colv.c_str()));
  }
  source.close();
  return v;
}

std::vector<int> read_my_tile_ids(int myRank,
				  const std::string & fileName)
{

  std::vector<int> v;
  std::ifstream source;
  source.open(fileName, std::ios_base::in);
  std::string line, tmp;

  int numCols = {};
  int count=0;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    if (count++ == 0){
      in >> tmp;
      numCols = std::atoi(tmp.c_str());
    }

    else{
      in >> tmp; const int rank= std::atoi(tmp.c_str());
      if (myRank==rank){
	for (int i=0; i<numCols-1; ++i){
	  in >> tmp; const int val = std::atoi(tmp.c_str());
	  if (val != -1){ v.push_back(val); }
	}
      }
    }
  }
  source.close();

  return v;
}

template<class CommType>
int calculate_total_num_of_tiles(CommType comm,
				 const std::vector<int> & myTileIds){
  auto maxTileId = (*std::max_element(myTileIds.cbegin(), myTileIds.cend()));
  int result = {};
  Teuchos::reduceAll(*comm, Teuchos::EReductionType::REDUCE_MAX, 1, &maxTileId, &result);
  return result+1; // because tileId starts at 0
}

std::map<int, std::vector<int>> read_my_state_vec_rows(const std::vector<int> & myTileIds)
{
  using v_t = std::vector<int>;
  using r_t = std::map<int, v_t>;
  r_t res;

  for (auto & it : myTileIds){
    const auto fileName = "state_vec_rows_wrt_full_mesh_p_" + std::to_string(it) + ".txt";
    res[it] = create_stdvector_and_fill_from_ascii(fileName);
  }
  return res;
}


class MyApp{
protected:
  using map_t		= Tpetra::Map<>;
  using go_t		= typename map_t::global_ordinal_type;
  using lo_t		= typename map_t::local_ordinal_type;
  using tcomm_t		= Teuchos::MpiComm<int>;
  using rcpcomm_t	= Teuchos::RCP<const tcomm_t>;
  using rcpmap_t	= Teuchos::RCP<const map_t>;

  template<typename T> using stdrcp = std::shared_ptr<T>;
  using crs_graph_type = Tpetra::CrsGraph<>;

public:
  /* these types exposed because need to be detected */
  using scalar_type	= double;
  using state_type	= Tpetra::Vector<>;
  using velocity_type	= Tpetra::Vector<>;

public:
  MyApp(rcpcomm_t comm) : comm_{comm} {
    setup();
  }
  ~MyApp() = default;

public:
  auto dataMapRcp(){ return dataMap_; };

  velocity_type createVelocity() const{
    velocity_type R(dataMap_);
    R.putScalar(0);
    return R;
  }

  void velocity(const state_type & /*u*/,
		const scalar_type  /*t*/,
		velocity_type & rhs) const{
    rhs.randomize();
  }

protected:
  void setup(){
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    wrappedCout_ = getFancyOStream(rcpFromRef(std::cout));
    myRank_   = comm_->getRank();
    totRanks_ = comm_->getSize();

    Kokkos::View<go_t*> myGidsV;
    const auto fileName= "tpetra_map_gids_" + std::to_string(myRank_) + ".txt";
    read_integers_from_ascii_into_view(fileName, myGidsV);
    dataMap_ = Teuchos::rcp(new map_t(-1, myGidsV, 0, comm_));

    // if (myRank_ == 0){
    //   std::cout << myRank_ << " "
    // 		<< totRanks_ << " "
    // 		<< dataMap_->getGlobalNumElements()
    // 		<< std::endl;
    // }

    NumMyElem_ = dataMap_->getNodeNumElements();
    dataMap_->describe(*wrappedCout_, Teuchos::VERB_LOW); //EXTREME);
  };

protected:
  Teuchos::RCP<Teuchos::FancyOStream> wrappedCout_;
  rcpcomm_t comm_{};
  rcpmap_t dataMap_{};
  int myRank_{};
  int totRanks_{};
  lo_t NumMyElem_{};
};

template<class TpetraMapType>
std::vector<int>
find_tile_bases_row_indices_that_I_own(const TpetraMapType & dataMap,
				       const std::vector<int> & tileStateVecRows)
{
  std::vector<int> result;
  auto myGlobalInd = dataMap.getMyGlobalIndices();
  for (int i=0; i<tileStateVecRows.size(); ++i){
    if (dataMap.isNodeGlobalElement( tileStateVecRows[i] )){
      result.push_back(i);
    }
  }

  return result;
}

struct ModesDataHolder
{
  using basis_mat_type = Kokkos::View<double**>;

  int rank_ = {};
  int modesPerTile_ = {};
  const std::vector<int> & myTileIds_;
  std::map<int, std::vector<int>> tileLocalIndices_;
  std::map<int, basis_mat_type> phi_;

  template<class TpetraMapType>
  ModesDataHolder(int rank, int modesPerTile,
		  const std::vector<int> & myTileIds,
		  const std::map<int, std::vector<int>> & myStateVecRowsForEachTile,
		  const TpetraMapType & dataMap)
    : rank_(rank), modesPerTile_(modesPerTile),
      myTileIds_(myTileIds)
  {
    auto myGlobalInd = dataMap.getMyGlobalIndices();

    for (int i=0; i<myTileIds_.size(); ++i){
      const int tileId = myTileIds[i];
      auto rows = find_tile_bases_row_indices_that_I_own(dataMap,
							 myStateVecRowsForEachTile.at(tileId));
      tileLocalIndices_[tileId] = rows;
      std::cout << "rank = " << " " << rank << " "
		<< "tileId = " << " " << tileId << " "
		<< "rows = " << rows.size() << '\n';

      basis_mat_type M("M", rows.size(), modesPerTile);
      const auto currFile = "bases_p_" + std::to_string(tileId) + ".txt";
      read_phi_from_ascii_into_view(currFile, M, rows);
      phi_[tileId] = M;
      // if (rank==1){
      // 	std::cout << tileId << " ";
      // 	for (int j=0; j<M.extent(1); ++j){ std::cout << M(0,j) << " "; }
      // 	for (int j=0; j<M.extent(1); ++j){ std::cout << M(M.extent(0)-1,j) << " "; }
      // 	std::cout << "\n";
      // }
      // if (rank==1 && tileId==6){
      // 	std::for_each(rows.begin(), rows.end(), [](int v){ std::cout << v << '\n'; });
      // }
    }

    // // figure out what is largest count of stencil dofs among tiles I own
    // int maxCount=0;
    // for (auto it : myStateVecRowsForEachTile){
    //   maxCount = std::max(maxCount, (int) it.second.size());
    // }
    // std::cout << "rank = " << " " << rank << " "
    // 	      << "maxCount = " << " " << maxCount << '\n';

    // // read bases
    // Kokkos::resize(phi_, myTileIds_.size(), maxCount, modesPerTile_);
    // for (int i=0; i<myTileIds_.size(); ++i){
    //   const int tileId = myTileIds[i];
    //   auto sv = Kokkos::subview(phi_, i, Kokkos::ALL(), Kokkos::ALL());
    //   const auto currFile = "bases_p_" + std::to_string(tileId) + ".txt";
    //   read_phi_from_ascii_into_view(currFile, sv);
    // }
  }

  const std::vector<int> & tileLocalIndices(int tileId) const {
    return tileLocalIndices_.at(tileId);
  }

  basis_mat_type viewBasis(int tileId) const{
    return phi_.at(tileId);
  }
};

struct OdRomDecoderJacobian
{
  using fom_state_type = Tpetra::Vector<>;

  OdRomDecoderJacobian(int rank, int modesPerTile,
		       const std::vector<int> & myTileIds,
		       const std::map<int, std::vector<int>> & myStateVecRowsForEachTile,
		       const ModesDataHolder & modesHolder)
    : rank_(rank), modesPerTile_(modesPerTile),
      myTileIds_(&myTileIds),
      myStateVecRowsForEachTile_(&myStateVecRowsForEachTile),
      modesHolder_(&modesHolder)
  {}

  template <class DataType, class ...Props>
  void applyTranspose(Tpetra::Vector<> x,
		      const Kokkos::View<DataType, Props...> y) const
  {
    // y = A^T x

    // here we need to be careful because after this operation
    // y should be same for all ranks.
    // but we handle partitions that exist only on some ranks
    // so to do this properly we need to first to a LOCAL operation
    // and store the result locally, which will yield a local
    // y with some zeros and then we need to do a global communication
    // to merge together all the y's from all the ranks

    using go_t		= typename decltype(x)::global_ordinal_type;
    using lo_t		= typename decltype(x)::local_ordinal_type;
    auto comm = x.getMap()->getComm();
    auto locallyRepMap = Tpetra::createLocalMap<lo_t, go_t>(y.extent(0), comm);
    Tpetra::Vector<> replY(locallyRepMap);

    auto replY_h0 = replY.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    auto replY_h  = Kokkos::subview(replY_h0, Kokkos::ALL(), 0);

    for (int tileIt=0; tileIt<myTileIds_->size(); ++tileIt){
      const int tileId = (*myTileIds_)[tileIt];

      auto phi = modesHolder_->viewBasis(tileId);

      const int romStateSpanStart = tileId*modesPerTile_;
      const int romStateSpanEnd   = romStateSpanStart + modesPerTile_;
      auto ysp = Kokkos::subview(replY_h, std::make_pair(romStateSpanStart, romStateSpanEnd));

      const auto & stateVecRows = myStateVecRowsForEachTile_->at(tileId);
      const auto & tileLocIndices = modesHolder_->tileLocalIndices(tileId);
      auto x_h = x.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
      const auto x_map = x.getMap();

      for (int i=0; i<tileLocIndices.size(); ++i){
	const int rowIndex = tileLocIndices[i];
	const auto gid = stateVecRows[rowIndex];
	const auto lid = x_map->getLocalElement(gid);
	for (int k=0; k<ysp.extent(0); ++k){
	  ysp(k) += phi(i, k) * x_h(lid,0);
	}
      }
    }

    replY.reduce();
    for (int k=0; k<y.extent(0); ++k){
      y(k) = replY_h(k);
    }
  }

private:
  int rank_ = {};
  int modesPerTile_ = {};
  const std::vector<int> * myTileIds_;
  const std::map<int, std::vector<int>> * myStateVecRowsForEachTile_;
  const ModesDataHolder * modesHolder_;
};

class OdRomDecoder
{
public:
  using fom_state_type = Tpetra::Vector<>;
  using jacobian_type  = OdRomDecoderJacobian;

  OdRomDecoder(int rank, int modesPerTile,
	       const std::vector<int> & myTileIds,
	       const std::map<int, std::vector<int>> & myStateVecRowsForEachTile,
	       const ModesDataHolder & modesHolder)
    : rank_(rank), modesPerTile_(modesPerTile),
      myTileIds_(myTileIds),
      myStateVecRowsForEachTile_(myStateVecRowsForEachTile),
      modesHolder_(modesHolder),
      jac_(rank, modesPerTile, myTileIds, myStateVecRowsForEachTile, modesHolder)
  {}

  template <class DataType, class ...Props>
  void applyMapping(const Kokkos::View<DataType, Props...> romState,
                    fom_state_type & fomState) const
  {

    for (int i=0; i<myTileIds_.size(); ++i){
      const int tileId = myTileIds_[i];

      const int romStateSpanStart = tileId*modesPerTile_;
      const auto spbounds = std::make_pair(romStateSpanStart,
					   romStateSpanStart + modesPerTile_);
      auto sp = Kokkos::subview(romState, std::move(spbounds));
      // if (rank_==1){
      // 	std::cout << tileId << " ";
      // 	for (int j=0; j<sp.extent(0); ++j){
      //          std::cout << sp(j) << " "; }
      // 	std::cout << "\n";
      // }

      auto phi = modesHolder_.viewBasis(tileId);
      // std::cout << rank_ << " " << tileId << " "
      // 	   << phi.extent(0) << " "
      //           << phi.extent(1) << "\n";

      Kokkos::View<double*> tmp("tmp", phi.extent(0));
      KokkosBlas::gemv("N", 1., phi, sp, 0., tmp);

      const auto & stateVecRows = myStateVecRowsForEachTile_.at(tileId);
      const auto & tileLocIndices = modesHolder_.tileLocalIndices(tileId);
      for (int k=0; k<tileLocIndices.size(); ++k){
	const int rowIndex = tileLocIndices[k];
	const auto gid = stateVecRows[rowIndex];
	fomState.replaceGlobalValue(gid, tmp(k));
      }
      // if (rank_==1 && tileId==2){
      // 	std::for_each(tileLocIndices.begin(),
      //                      tileLocIndices.end(),
      //                      [](int v){ std::cout << v << '\n'; });
      // }
    }
  }

  // if applicable, update the Jacobian for a given state
  template <typename OperandType>
  void updateJacobian(const OperandType & romOperand) {}

  // return a const reference to the Jacobian object
  const jacobian_type & jacobianCRef() const{ return jac_; }

private:
  int rank_ = {};
  int modesPerTile_ = {};
  const std::vector<int> & myTileIds_;
  const std::map<int, std::vector<int>> & myStateVecRowsForEachTile_;
  const ModesDataHolder & modesHolder_;
  OdRomDecoderJacobian jac_;
};

#include "pressio/ops.hpp"

namespace pressio{ namespace ops{

template<class ScalarType, class DataType, class ...Props>
void product(::pressio::transpose /*mode*/, ScalarType alpha,
	     const OdRomDecoderJacobian & A,
	     Tpetra::Vector<> x, ScalarType beta,
	     const Kokkos::View<DataType, Props...> y)
{
  // y = beta * y + alpha * A^T x

  using pressio::utils::Constants;
  if (alpha == Constants<ScalarType>::one() &&
      beta  == Constants<ScalarType>::zero()){
    A.applyTranspose(x, y);
  }
  else{
    throw std::runtime_error("missing impl");
  }
}

}} //end namespace pressio::ops

#include "pressio/rom_galerkin.hpp"

TEST(odrom, tpetra_draft)
{
  using tcomm_t	  = Teuchos::MpiComm<int>;
  using rcpcomm_t = Teuchos::RCP<const tcomm_t>;

  // namespace plog = pressio::log;
  // plog::initialize(pressio::logto::terminal);
  {

    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    rcpcomm_t comm = Teuchos::rcp (new tcomm_t(MPI_COMM_WORLD));

    MyApp appObj(comm);
    auto dataMapRcp = appObj.dataMapRcp();
    std::cout << std::flush;

    const auto myTileIds = read_my_tile_ids(myRank, "rank_partition_mapping.txt");
    sleep((myRank+1.)*0.5);
    std::cout << "rank: " << myRank << " tileIds: ";
    std::for_each(myTileIds.begin(), myTileIds.end(), [](int v){ std::cout << v << ' '; });
    std::cout << "\n";

    const auto myStateVecRowsForEachTile = read_my_state_vec_rows(myTileIds);
    const int totalNumTiles = calculate_total_num_of_tiles(comm, myTileIds);
    std::cout << "rank: " << totalNumTiles << "\n";
    std::cout << std::flush;

    int K = 3;
    ModesDataHolder modesHolder(myRank, K, myTileIds,
				myStateVecRowsForEachTile, *dataMapRcp);

    OdRomDecoder decoder(myRank, K, myTileIds,
			 myStateVecRowsForEachTile, modesHolder);

    using rom_type = Kokkos::View<double*>;
    rom_type romState("romState", K*totalNumTiles);
    for (int iTile=0; iTile<totalNumTiles; ++iTile){
      const int start = iTile*K;
      auto sp = Kokkos::subview(romState, std::make_pair(start, start+K));
      for (int j=0; j<sp.extent(0); ++j){ sp(j) = (double) (iTile+1); }
    }
    if (myRank==0){
      write_kokkos_view_to_ascii_file("rom_state_initial.txt", romState);
    }


    // reconstruction
    typename MyApp::state_type fomTmpState(dataMapRcp);
    fomTmpState.putScalar(0.);
    decoder.applyMapping(romState, fomTmpState);
    write_tpetra_vector_to_ascii_file("fomTmpState_"+std::to_string(myRank)+".txt", fomTmpState);
    std::cout << std::flush;

    {
      // projection
      auto velo = appObj.createVelocity();
      velo.putScalar(0.);
      appObj.velocity(fomTmpState, 0.0, velo);
      write_tpetra_vector_to_ascii_file("velo_"+std::to_string(myRank)+".txt", velo);
      const auto & dJ = decoder.jacobianCRef();
      pressio::ops::product(pressio::transpose(), 1., dJ, velo, 0., romState);
      //if (myRank==0){
	write_kokkos_view_to_ascii_file("rom_state_projected_"+std::to_string(myRank)+".txt", romState);
	//}
    }

    // typename MyApp::state_type fomRefState(dataMapRcp);
    // namespace pode = pressio::ode;
    // namespace pgal = pressio::rom::galerkin;
    // constexpr auto odeScheme = pode::StepScheme::RungeKutta4;
    // auto galProb = pgal::create_default_explicit_problem(odeScheme, appObj, decoder,
    // 							 romState, fomRefState);
    // pode::advance_n_steps(galProb, romState, 0., 0.01, 10);

    std::cout << std::flush;
    comm->barrier();
  }
  sleep(3.);
  //plog::finalize();
}
