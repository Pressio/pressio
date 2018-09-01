
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( TRILINOS
  # REQUIRED_HEADERS  SimpleTpl.hpp
  REQUIRED_LIBS_NAMES libamesos2.dylib libanasazi.dylib libanasaziepetra.dylib libanasazitpetra.dylib 
  libaztecoo.dylib libbelos.dylib libbelosepetra.dylib libbelostpetra.dylib libepetra.dylib 
  libepetraext.dylib libgaleri-epetra.dylib libifpack.dylib libifpack2.dylib libkokkosalgorithms.dylib 
  libkokkoscontainers.dylib libkokkoscore.dylib libkokkoskernels.dylib libkokkostsqr.dylib 
  libloca.dylib liblocalapack.dylib libml.dylib libModeLaplace.dylib libnox.dylib 
  libnoxlapack.dylib librol.dylib libteuchoscomm.dylib libteuchoscore.dylib libteuchoskokkoscomm.dylib 
  libteuchoskokkoscompat.dylib libteuchosnumerics.dylib libteuchosparameterlist.dylib 
  libteuchosparser.dylib libteuchosremainder.dylib libtpetra.dylib libtpetraclassic.dylib 
  libtpetraclassiclinalg.dylib libtpetraclassicnodeapi.dylib libtpetraext.dylib 
  libtpetrainout.dylib libtrilinosss.dylib libtriutils.dylib

  MUST_FIND_ALL_HEADERS
  )




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
