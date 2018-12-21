#ifndef SOLVERS_LINEAR_TAGS_HPP
#define SOLVERS_LINEAR_TAGS_HPP

#include "solvers_ConfigDefs.hpp"

namespace rompp{ namespace solvers{ namespace linear {

// // Linear dense solvers types
// struct ColPivHouseholderQR {};
// struct CompleteOrthogonalDecomposition {};

// Linear iterative solvers types
struct CG {};
struct LSCG {};
struct Bicgstab {};
// struct Gmres {};

// // Preconditioner types
// struct Jacobi {};
// struct DefaultPreconditioner {};

}}}//end namespace rompp::solvers::linear



#endif
