
#ifdef HAVE_TRILINOS
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_POLICIES_TRILINOS_HPP_
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_POLICIES_TRILINOS_HPP_

#include "AztecOO.h"

namespace pressio{
namespace solvers{
namespace trilinos_policies {


// Linear solvers policies
struct CG {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_cg);
  }
};


struct Gmres {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_gmres);
  }
};


struct BicgStab {
  static void set_solver_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_solver, AZ_bicgstab);
  }
};


// Preconditioners policies
struct DefaultPreconditioner {
  static void set_preconditioner_type(AztecOO& solver) {
    return;
  }
};

struct DecompIlu {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
  }
};


struct DecompIlut {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver.SetAztecOption(AZ_overlap, 1);
    solver.SetAztecOption(AZ_ilut_fill, 3.0);
  }
};


struct DecompIcc {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_icc); 
  }
};


struct Jacobi {
  static void set_preconditioner_type(AztecOO& solver) {
    solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    solver.SetAztecOption(AZ_poly_ord, 1);
  }
};


} // end namespace trilinos_policies
} // end namespace solvers
}//end namespace pressio
#endif
#endif
