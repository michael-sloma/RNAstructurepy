##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'Use "make TurboFold" to create executable "TurboFold."'
	@echo 'Use "make TurboFold_Library" to make the TurboFold library.'
	@echo 'Use "make TurboFold-smp" to create executable "TurboFold-smp."'
	@echo 'Use "make TurboFoldSMP_Library" to make the TurboFold SMP library.'

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=..
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

##########
## Define targets.
##########

# Build the TurboFold serial text interface.
TurboFold: ../exe/TurboFold
../exe/TurboFold: ${ROOTPATH}/TurboFold/TurboFold_Interface.o ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${TURBOFOLD_SERIAL_FILES}
	${LINK} ${ROOTPATH}/TurboFold/TurboFold_Interface.o ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${TURBOFOLD_SERIAL_FILES}

# Build the TurboFold shared library for serial calculations.
TurboFold_Library: ${TURBOFOLD_SERIAL_FILES}
	ar rc ${ROOTPATH}/exe/${TURBOFOLD_LIBRARY} ${TURBOFOLD_SERIAL_FILES}

# Build the TurboFold SMP text interface.
TurboFold-smp: ../exe/TurboFold-smp
../exe/TurboFold-smp: ${ROOTPATH}/TurboFold/TurboFold_Interface-smp.o ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${TURBOFOLD_SMP_FILES}
	${LINK_IN_EXE_DIR} ${ROOTPATH}/TurboFold/TurboFold_Interface-smp.o ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${TURBOFOLD_SMP_FILES} -lpthread

# Build the TurboFold shared library for SMP calculations.
TurboFoldSMP_Library: ${TURBOFOLD_SMP_FILES}
	ar rc ${ROOTPATH}/exe/${TURBOFOLD_SMP_LIBRARY} ${TURBOFOLD_SMP_FILES}

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

clean:
	rm -f *.o
	rm -f ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER}
	rm -f ${TURBOFOLD_SERIAL_FILES}
	rm -f ${TURBOFOLD_SMP_FILES}

realclean: clean
	rm -f ${ROOTPATH}/exe/TurboFold
	rm -f ${ROOTPATH}/exe/${TURBOFOLD_LIBRARY}
	rm -f ${ROOTPATH}/exe/TurboFold-smp
	rm -f ${ROOTPATH}/exe/${TURBOFOLD_SMP_LIBRARY}

