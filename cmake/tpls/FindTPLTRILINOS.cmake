
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


  



# libkokkosalgorithms.dylib libkokkoscontainers.dylib 
# libkokkoscore.dylib libkokkoskernels.dylib libkokkostsqr.dylib 
# libamesos2.dylib libanasazi.dylib libanasaziepetra.dylib 
# libanasazitpetra.dylib libaztecoo.dylib libbelos.dylib libbelosepetra.dylib 
# libbelostpetra.dylib libepetra.dylib libepetraext.dylib 
# libgaleri-epetra.dylib libifpack.dylib libifpack2.dylib 
# libloca.dylib liblocalapack.dylib libml.dylib 
# libModeLaplace.dylib libnox.dylib libnoxlapack.dylib 
# librol.dylib libteuchoscomm.dylib libteuchoscore.dylib 
# libteuchoskokkoscomm.dylib libteuchoskokkoscompat.dylib 
# libteuchosnumerics.dylib libteuchosparameterlist.dylib 
# libteuchosparser.dylib libteuchosremainder.dylib 
# libtpetra.dylib libtpetraclassic.dylib libtpetraclassiclinalg.dylib 
# libtpetraclassicnodeapi.dylib libtpetraext.dylib 
# libtpetrainout.dylib libtrilinosss.dylib libtriutils.dylib



# INCLUDE(TribitsTplDeclareLibraries)

# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Trilinos
#   REQUIRED_LIBS_NAMES xpetra-ext tpi tpetrainout tpetraext
# tpetraclassiclinalg tpetraclassic thyraepetraext
# teuchoskokkoscomm stratimikosbelos shards
# sacado rol pamgen_extras pamgen muelu-adapters
# locathyra localapack locaepetra loca kokkostsqr
# kokkoscontainers kokkosalgorithms intrepid ifpack2
# ifpack2-adapters gtest galeri-xpetra belostpetra
# belosepetra anasazitpetra anasaziepetra anasazi
# amesos2 ModeLaplace zoltan2 xpetra-sup thyratpetra
# noxlapack noxepetra nox muelu-interface teko
# stratimikos muelu isorropia belos xpetra
# tpetrakernels tpetra stratimikosml stratimikosifpack
# stratimikosaztecoo stratimikosamesos ml ifpack
# aztecoo amesos zoltan tpetraclassicnodeapi thyraepetra
# thyracore teuchosremainder teuchosnumerics
# teuchoskokkoscompat rtop epetraext triutils galeri 
# epetra teuchoscomm teuchosparameterlist teuchoscore kokkoscore
# )
