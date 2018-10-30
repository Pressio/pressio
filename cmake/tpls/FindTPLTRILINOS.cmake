
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( TRILINOS
  REQUIRED_LIBS_NAMES 
  libkokkoskernels.dylib #kokkos kernels
  libkokkosalgorithms.dylib libkokkoscontainers.dylib #kokkos
  libteuchoskokkoscomm.dylib libteuchoskokkoscompat.dylib #kokkos
  libkokkoscore.dylib  libkokkostsqr.dylib #kokkos
  libteuchoscomm.dylib libteuchoscore.dylib libteuchosnumerics.dylib #teuchos 
  libteuchosparameterlist.dylib libteuchosparser.dylib libteuchosremainder.dylib #teuchos
  libepetra.dylib libgaleri-epetra.dylib #epetra 
  libtpetra.dylib libtpetraclassic.dylib
  libtpetraclassiclinalg.dylib libtpetraclassicnodeapi.dylib #tpetra
  libtpetraext.dylib libtpetrainout.dylib #tpetra
  libtriutils.dylib
  libtrilinosss.dylib #trilinosSS
  libepetraext.dylib #epetra ext
  libxpetra-sup.dylib libxpetra.dylib #xpetra
  libaztecoo.dylib #aztecoo  
  libgaleri-xpetra.dylib
  libifpack.dylib
  libml.dylib
  libbelos.dylib libbelosepetra.dylib libbelostpetra.dylib libbelosxpetra.dylib #belos
  libamesos2.dylib
  libanasazi.dylib libanasaziepetra.dylib libanasazitpetra.dylib
  libifpack2.dylib
  libnox.dylib libnoxepetra.dylib libnoxlapack.dylib
  libmuelu-adapters.dylib libmuelu-interface.dylib  libmuelu.dylib
  librol.dylib
  MUST_FIND_ALL_HEADERS
  )




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
