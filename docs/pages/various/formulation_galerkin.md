
# The Galerkin ROM formulation

Galerkin projection can be derived by 
minimizing the *time-continuous* residual over the trial manifold.
The resulting model can be obtained by substituting the approximate state
and the corresponding velocity @f$f(\tilde{x})@f$  into the time-continuous residual and minimizing 
its (weighted) @f$\ell^2@f$-norm,
which yields a sequence of
residual-minimization problems
@m_class{m-default}

@f[
	\dot{\hat{x}}(t;\mu)  =
	\underset{\xi \in R^{p}}{arg min}
    \left\|
	A \left(J(\hat{x}(t;\mu))\xi - f\left(x_{ref}(\mu)+g(\hat{x}(t;\mu);\mu\right) \right)
    \right\|_2^2
@f]

which can be written equivalently as

@f[
	\dot{\hat{x}}(t;\mu)  = \left(AJ(\hat{x}(t;\mu) \right)^+ f\left(x_{ref}(\mu)+g(\hat{x}(t;\mu);\mu\right)
@f]

where the superscript + denotes the Moore-Pentrose pseudoinverse and @f$\hat{x}(0;\mu)=\hat{x}^0(\mu)@f$ is the reduced initial condition.
The matrix @f$A \in R^{z \times N}@f$
with @f$z \leq N@f$ denotes a weighting
matrix that can enable hyper-reduction.




