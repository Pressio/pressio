
#ifndef PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_MASKERS_HPP_
#define PRESSIO_TEST_ROM_GALERKIN_MASKED_CORRECT_COMMON_MASKERS_HPP_

#include <gtest/gtest.h>

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "pressio/ops.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#endif

// for explicit, masker acts on FOM velicity
struct MaskerExplicitEigen
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = typename TrivialFomOnlyVelocityEigen::velocity_type;

  MaskerExplicitEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(sample_indices_.size());
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

// for implicit, masker acts on FOM velicity and FOM apply jac result
struct MaskerImplicitEigen
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;

  MaskerImplicitEigen(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.cols());
  }

  template<class TimeType>
  void operator()(const vec_operand_type & operand, TimeType time, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  template<class TimeType>
  void operator()(const mat_operand_type & operand, TimeType time, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.cols(); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

// for explicit, masker acts on FOM velicity
template<class ScalarType>
struct MaskerExplicitCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = ::pressiotests::MyCustomVector<ScalarType>;

  MaskerExplicitCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(sample_indices_.size());
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};

// for implicit, masker acts on FOM velicity and FOM apply jac result
template<class ScalarType>
struct MaskerImplicitCustomTypes
{
  const std::vector<int> sample_indices_ = {};
  using vec_operand_type = ::pressiotests::MyCustomVector<ScalarType>;
  using mat_operand_type = ::pressiotests::MyCustomMatrix<ScalarType>;

  MaskerImplicitCustomTypes(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  vec_operand_type createApplyMaskResult(const vec_operand_type & operand) const{
    return vec_operand_type(sample_indices_.size());
  }

  mat_operand_type createApplyMaskResult(const mat_operand_type & operand) const{
    return mat_operand_type(sample_indices_.size(), operand.extent(1));
  }

  template<class TimeType>
  void operator()(const vec_operand_type & operand, TimeType time, vec_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }

  template<class TimeType>
  void operator()(const mat_operand_type & operand, TimeType time, mat_operand_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (std::size_t j=0; j<operand.extent(1); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// for explicit, masker acts on FOM velicity
struct MaskerExplicitTpetra
{
  using tcomm = Teuchos::Comm<int>;
  using map_t = Tpetra::Map<>;
  using GO = typename map_t::global_ordinal_type;
  using LO = typename map_t::local_ordinal_type;

  Teuchos::RCP<const tcomm> comm_;
  Teuchos::RCP<const map_t> maskedMap_;
  std::vector<int> sample_indices_ = {};
  using operand_type = typename TrivialFomOnlyVelocityTpetra::velocity_type;

  MaskerExplicitTpetra(Teuchos::RCP<const tcomm> comm, 
                       Teuchos::RCP<const map_t> full_map,                       
                       std::vector<int> sample_indices) 
  : comm_(comm), sample_indices_(sample_indices)
  {
    auto mygids_fm = full_map->getMyGlobalIndices();
    // 100 is large enough for this test, only a subset 
    // will be picked by the map constructor
    GO mymaskedGids[100];
    GO count=0;
    for (std::size_t i=0; i<mygids_fm.extent(0); ++i)
    {
      const bool r = std::find(sample_indices_.cbegin(), sample_indices_.cend(), mygids_fm(i))!=sample_indices_.cend();
      if (r){ mymaskedGids[count++] = mygids_fm(i); }
    }

    maskedMap_ = Teuchos::rcp(new map_t(-1, mymaskedGids, count, 0, comm_));
  }

  Teuchos::RCP<const map_t> map(){ return maskedMap_; }

  operand_type createApplyMaskResult(const operand_type & operand) const{
    return operand_type(maskedMap_);
  }

  template<class TimeType>
  void operator()(const operand_type & operand, TimeType time, operand_type & result) const
  {
    using importer_t = Tpetra::Import<LO, GO>;
    importer_t importer(operand.getMap(), maskedMap_);
    result.doImport(operand, importer, Tpetra::CombineMode::REPLACE);
  }
};
#endif

#endif
