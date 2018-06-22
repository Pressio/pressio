TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  core   		packages/core  		PT
  solvers 		packages/solvers 	PT
  ode    		packages/ode   		PT
  apps   		packages/apps  		PT
  #
  # optimization	packages/optimization   PT
  #svd    	packages/svd   PT
  #rom    	packages/rom   PT
  )


# potential keywords: PT, ST

#TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(WrapExternal Windows)
#TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(InsertedPkg)

