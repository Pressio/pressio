
#ifndef ROM_IMPL_GALERKIN_HELPERS_HPP_
#define ROM_IMPL_GALERKIN_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class T = void>
void valid_scheme_for_explicit_galerkin_else_throw(::pressio::ode::StepScheme name,
						   const std::string & str){
  if (!::pressio::ode::is_explicit_scheme(name)){
    throw std::runtime_error(str + " requires an explicit stepper");
  }
}

template<class T = void>
void valid_scheme_for_implicit_galerkin_else_throw(::pressio::ode::StepScheme name,
						   const std::string & str){
  if (!::pressio::ode::is_implicit_scheme(name)){
    throw std::runtime_error(str + " requires an implicit stepper");
  }
}

// --------------------------------------------------------------
// CreateGalerkinRhs
// --------------------------------------------------------------
template<class T, class = void> struct CreateGalerkinRhs;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct CreateGalerkinRhs<
  T, mpl::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  >
{
  T operator()(std::size_t ext){ return T(ext); }
};
#endif

// --------------------------------------------------------------
// CreateGalerkinMassMatrix
// --------------------------------------------------------------
template<class T, class = void>
struct CreateGalerkinMassMatrix;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct CreateGalerkinMassMatrix<
  T, mpl::enable_if_t< ::pressio::is_dense_matrix_eigen<T>::value >
  >
{
  T operator()(std::size_t ext){ return T(ext, ext); }
};
#endif

// ------------------------------------------
// this is an alias becuase it can done in the same way
template<class JacType, class = void>
using CreateGalerkinJacobian = CreateGalerkinMassMatrix<JacType>;

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_GALERKIN_HELPERS_HPP_
