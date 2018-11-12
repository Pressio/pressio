
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( TRILINOS
  REQUIRED_LIBS_NAMES
  kokkoskernels 
  kokkosalgorithms
  kokkoscontainers 
  teuchoskokkoscomm
  teuchoskokkoscompat 
  kokkoscore
  kokkostsqr 
  teuchoscomm
  teuchoscore
  teuchosnumerics 
  teuchosparameterlist
  teuchosparser
  teuchosremainder 
  epetra
  tpetra
  tpetraclassic
  tpetraclassiclinalg
  tpetraclassicnodeapi 
  tpetraext
  tpetrainout
  epetraext
  aztecoo 
  anasazi
  anasaziepetra
  anasazitpetra
  ifpack
  ifpack2
  galeri-epetra
  
  MUST_FIND_ALL_HEADERS
  )



  #libgaleri-xpetra.dylib
  #libml.dylib
  #libbelos.dylib libbelosepetra.dylib libbelostpetra.dylib libbelosxpetra.dylib #belos
  #libamesos2.dylib
  #libnox.dylib libnoxepetra.dylib libnoxlapack.dylib
  #libmuelu-adapters.dylib libmuelu-interface.dylib  libmuelu.dylib
  #librol.dylib
