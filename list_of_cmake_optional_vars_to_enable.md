
This file contains a list of cmake options (and their description)
that can be turned on/off by the user to enable/disable certain features
at compile time since these are used for pre-processor directives.

- PRESSIO_ENABLE_DEBUG_PRINT
	- description:	turn on/off print statemetns for debugging
	- enable with: 	-DPRESSIO_ENABLE_DEBUG_PRINT=[ON/OFF]
	- prerequisites: 	none
	- default: 		OFF


- PRESSIO_ENABLE_TPL_MPI
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_MPI=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


- PRESSIO_ENABLE_TPL_TRILINOS
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_TRILINOS=[ON/OFF]
	- prerequisites: 	need to have Trilinos installed if you want this on
	- default: OFF


- PRESSIO_ENABLE_TPL_KOKKOS
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_KOKKOS=[ON/OFF]
	- prerequisites: 	need to have Trilinos or Kokkos installed if you want this on
	- default: 		OFF
	- It is set ON automatically if -PRESSIO_ENABLE_TPL_TRILINOS=ON


- PRESSIO_ENABLE_TEUCHOS_TIMERS
	- description: 	turn off/on teuchos timers
	- prerequisites: 	DPRESSIO_ENABLE_TPL_TRILINOS=ON
	- ON if -DPRESSIO_ENABLE_TPL_TRILINOS=ON and PRESSIO_ENABLE_DEBUG_PRINT=ON

- PRESSIO_ENABLE_TPL_BLAS
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_BLAS=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


- PRESSIO_ENABLE_TPL_LAPACK
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_LAPACK=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


- PRESSIO_ENABLE_TPL_PYBIND11
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_PYBIND11=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


- PRESSIO_ENABLE_TPL_BLAZE
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_BLAZE=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


- PRESSIO_ENABLE_TPL_ARMADILLO
	- description: 	self-explanatory
	- enable with: 	-DPRESSIO_ENABLE_TPL_ARMADILLO=[ON/OFF]
	- prerequisites: 	need to have it installed if you want this on
	- default: 		OFF


