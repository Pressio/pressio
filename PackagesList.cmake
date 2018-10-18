

TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  # make sure the order in which you list packages here makes
  # sense in terms of dependencies. packages listed lower should
  # depend on those lister higher. Dependencies should be read
  # from top to bottom. 
  core packages/core PT
  solvers packages/solvers PT
  qr packages/qr PT
  svd packages/svd PT
  ode packages/ode PT
  optimization packages/optimization PT
  rom packages/rom PT
  )


# potential keywords: PT, ST
#TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
#TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(InsertedPkg)

