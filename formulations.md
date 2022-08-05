
# STEADY GALERKIN

## Default

- equation: phi^T f(phi x) = 0
- R = phi^T f(phi x)
- J = phi^T df/dx(phi x) phi
- should always accept trialSpace/affineTrialSpace

## Hypred

- equation: MJOP f(phi x) = 0
- R = MJOP f(phi x)
- J = MJOP df/dx(phi x) phi
- should always accept trialSpace/affineTrialSpace

## Masked

- MJOP masked[f(phi x)] = 0
- R = MJOP masked[f(phi x)]
- J = MJOP masked[df/dx(phi x) phi]
- should always accept trialSpace/affineTrialSpace


<!-- ============================================================== -->


# UNSTEADY GALERKIN

## Default

- should always accept trialSpace/affineTrialSpace

- "phi^T M phi dx/dt = phi^T f"
  - M is time dependent
  - or M constant
- "dx/dt = phi^T f"

- mass matrix or not is inferred from the system API


# Hypred

- should always accept trialSpace/affineTrialSpace
- only for M constant or no M
- "phi^T M phi dx/dt = MJOP f"
- "dx/dt = MJOP f"

# Masked

- should always accept trialSpace/affineTrialSpace
- figure it out


<!-- ============================================================== -->


# STEADY LSPG

## Default

- arg min || r ||

## Hypred

- arg min ||A r_sample||
- A is used in the solver NOT in the actual LSPG class

## Masked

- arg min ||A r_sample||
- A is used in the solver NOT in the actual LSPG class


<!-- ============================================================== -->


# UNSTEADY LSPG

## Default

- arg min ||r||

## Hypred

- arg min ||A r_sample||
- A is used in the solver NOT in the actual LSPG class

## Masked

- arg min ||A r_sample||
- A is used in the solver NOT in the actual LSPG class
