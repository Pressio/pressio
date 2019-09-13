
TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  # TRIBITS wants the order in which you list packages here matching dependencies.
  # Dependencies should be read from top to bottom.
  # Packages listed lower should depend on those lister higher.
  mpl		packages/mpl PT
  utils		packages/utils PT
  containers	packages/containers PT
  apps		packages/apps PT
  qr		packages/qr PT
  svd		packages/svd PT
  optimization	packages/optimization PT
  solvers	packages/solvers PT
  ode		packages/ode PT
  rom		packages/rom PT
  )


# potential keywords: PT, ST
#TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
#TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(InsertedPkg)
