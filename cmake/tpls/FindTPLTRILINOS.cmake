
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( TRILINOS
  REQUIRED_LIBS_NAMES
  kokkosalgorithms
  kokkoscontainers
  kokkoscore
  teuchoskokkoscomm
  teuchoskokkoscompat
  teuchosremainder
  teuchosnumerics
  teuchoscomm
  teuchosparameterlist
  teuchosparser
  teuchoscore
  epetra
  epetraext
  ifpack
  aztecoo
  tpetraext
  tpetrainout
  tpetra
  kokkostsqr
  tpetraclassiclinalg
  tpetraclassicnodeapi
  tpetraclassic
  kokkoskernels
  ifpack2
  #anasazitpetra
  #anasaziepetra
  #anasazi
  kokkosalgorithms
  teuchosparameterlist

  MUST_FIND_ALL_HEADERS
  )
