
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
  anasazitpetra
  anasaziepetra
  anasazi
  kokkosalgorithms
  kokkoscontainers
  kokkoscore	
  teuchosparameterlist

  MUST_FIND_ALL_HEADERS
  )



  #libgaleri-xpetra.dylib
  #libml.dylib
  #libbelos.dylib libbelosepetra.dylib libbelostpetra.dylib libbelosxpetra.dylib #belos
  #libamesos2.dylib
  #libnox.dylib libnoxepetra.dylib libnoxlapack.dylib
  #libmuelu-adapters.dylib libmuelu-interface.dylib  libmuelu.dylib
  #librol.dylib
