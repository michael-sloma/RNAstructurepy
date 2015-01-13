##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'Use "make GUI" to create the RNAstructure Java interface.'

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=..
JAVA=yes
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

##########
## Define targets.
##########

# Rule to build the Java interface, C++ and Java code.
GUI: buildC buildJava

# Rule to build the C++ part of the RNAstructure Java interface only.
buildC: ${JAVA_LIBRARY_FILES} ${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator.o ${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator_wrap.o
	${LINK_PROXY_FOR_JAVA}${RNASTRUCTURE_LIBRARY} ${JAVA_LIBRARY_FILES} ${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator.o ${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator_wrap.o

# Rule to build the Java part of the RNAstructure Java interface only.
buildJava:
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_utilities/*.java
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_drawing/source/*/*.java
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_interface/source/*/*.java
	cd ${ROOTPATH}; jar cvmf RNAstructure_java_interface/Manifest.txt exe/RNAstructure.jar RNAstructure_java_interface/images/*.gif RNAstructure_java_interface/sounds/*.wav RNAstructure_java_interface/source/*/*.class RNAstructure_java_drawing/source/*/*.class RNAstructure_java_utilities/*.class

# Rule to build SWIG files.
# This rule is almost never necessary for the end user.
swigBuild:
	swig -java -c++ -package RNAstructure_java_interface.source.startup -outdir ${ROOTPATH}/RNAstructure_java_interface/source/startup ${ROOTPATH}/SWIG/RNAstructureBackendCalculator.i
	@rm -f ${ROOTPATH}/RNAstructure_java_interface/source/startup/*Proxy.java

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

# Rule to clean all possible objects.
clean:
	rm -f ${ROOTPATH}/RNAstructure_java_interface/source/*/*.class
	rm -f ${ROOTPATH}/RNAstructure_java_drawing/source/*/*.class
	rm -f ${ROOTPATH}/RNAstructure_java_utilities/*.class
	rm -f ${ROOTPATH}/RNAstructure_java_interface/SWIG/*.o
	rm -f ${ROOTPATH}/RNAstructure_java_drawing/SWIG/*.o
	rm -f ${HYBRID_FILES} ${DYNALIGN_SERIAL_FILES} ${OLIGO_FILES} ${TURBOFOLD_SERIAL_FILES} ${PLOT_HANDLER} ${STRUCTURE_IMAGER} ${ROOTPATH}/RNA_class/Multilign_object.o

# Rule to clean all possible objects and executables.
realclean: clean
	rm -f ${ROOTPATH}/exe/${RNASTRUCTURE_LIBRARY} ${ROOTPATH}/exe/RNAstructure.jar

