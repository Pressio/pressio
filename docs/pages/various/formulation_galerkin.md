
# The Galerkin ROM formulation

The Galerkin ROM is one of the most important techniques
for model reduction.

\todo put some references




@m_class{m-block m-default}

@parblock FOM FORMULATION

To formulate the problem, consider a dynamical system of the form
@f[
\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},t;\mathbf{\mu}), \qquad \mathbf{x}(0;\mu) = \mathbf{x}^0(\mu),
@f]
and let @f$\mathbf{r} \left(\dot{\mathbf{x}}, \mathbf{x}, t; \mathbf{\mu} \right) = \mathbf{0}@f$
be the (time-continuous) residual @f$\mathbf{r} \in R^{N}@f$
defined as @f$\mathbf{r} := \dot{\mathbf{x}}- \mathbf{f}(\mathbf{x},t;\mathbf{\mu})@f$.
@endparblock


\todo missing the approximation of the state


@m_class{m-block m-success}

@parblock
Galerkin projection can be derived by
minimizing the *time-continuous* residual over the trial manifold.
The resulting model can be obtained by substituting the approximate state
and the corresponding velocity @f$f(\tilde{x})@f$ into the time-continuous
residual and minimizing its (weighted) @f$\ell^2@f$-norm,
which yields a sequence of residual-minimization problems


@m_class{m-default}

@f[
	\dot{\hat{\mathbf{x}}}(t; \mathbf{\mu})  =
	\underset{\mathbf{\xi} \in R^{p}}{arg min}
    \left\|
	\mathbf{A} \left( \mathbf{J}(\hat{\mathbf{x}}(t;\mathbf{\mu}))\mathbf{\xi}
	- \mathbf{f}\left(\mathbf{x}_{ref}(\mathbf{\mu})
	+ \mathbf{g}(\hat{\mathbf{x}}(t;\mathbf{\mu});\mathbf{\mu}
	\right) \right)
    \right\|_2^2
@f]

which can be equivalently written as

@m_class{m-success}

@f[
	\dot{\hat{\mathbf{x}}}(t;\mathbf{\mu}) =
	\Big( \mathbf{A} \mathbf{J}(\hat{\mathbf{x}}(t;\mathbf{\mu}) \Big)^+
	\mathbf{A} \mathbf{f}
	\Big(\mathbf{x}_{ref}(\mathbf{\mu})
	+ \mathbf{g}(\hat{\mathbf{x}}(t;\mathbf{\mu}); \mathbf{\mu} \Big)
@f]

where the superscript + denotes the Moore-Pentrose pseudoinverse
and @f$\hat{\mathbf{x}}(0;\mathbf{\mu})=\hat{\mathbf{x}}^0(\mathbf{\mu})@f$
is the reduced initial condition. The matrix @f$\mathbf{A} \in R^{z \times N}@f$
with @f$z \leq N@f$ denotes a weighting matrix that can enable hyper-reduction.
@endparblock

\todo describe what kind of problems Galerkin is good for and those it is not good for


\todo put figure that shows how big the operators are


\todo Read [this](md_pages_adapter_apis_adapter_galerkin_api.html) to know more about what adapter API you need
to access the Galerkin ROM classes in pressio.
