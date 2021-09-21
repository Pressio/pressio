
# The LSPG ROM formulation

LSPG projection corresponds to minimizing
the (weighted) @f$\ell^2@f$-norm of the *time-discrete*
residual over the trial manifold.
Hence, the starting point for this approach is the residual
ODE formulation discretized in time with an arbitrary time-discretization method.

LSPG projection is derived by substituting the approximate state
in the time-discrete residual and minimizing its (weighted) @f$\ell^2@f$-norm,
which yields a sequence of
residual-minimization problems
@m_class{m-default}

@f[
	\hat{x}^n(\mu)  =
	\underset{\xi \in R^{p}}{arg min}
    \left\|
	A r^{n}\left(x_{ref}(\mu)+g(\xi);\mu)\right)
    \right\|_2^2,\quad
	n=1,\ldots,N_t
@f]

with initial condition @f$\hat{x}(0;\mu)=\hat{x}^0(\mu)@f$.
The matrix @f$A \in R^{z \times N}@f$
with @f$z \leq N@f$ denotes a weighting
matrix that can enable hyper-reduction.
