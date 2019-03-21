

TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  # make sure the order in which you list packages here makes
  # sense in terms of dependencies.
  # Dependencies should be read from top to bottom.
  # Packages listed lower should depend on those lister higher.
  core packages/core PT
  qr packages/qr PT
  solvers packages/solvers PT
  svd packages/svd PT
  ode packages/ode PT
  optimization packages/optimization PT
  rom packages/rom PT
  apps packages/apps PT
  )


# potential keywords: PT, ST
#TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
#TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(InsertedPkg)
