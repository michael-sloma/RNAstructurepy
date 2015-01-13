##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'Use "make plots" to create a self-contained dot plot drawing program.'
	@echo 'Use "make structures" to create a self-contained structure drawing program.'

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

# Rule to build the plot drawing program, C++ and Java code.
# This rule builds a dynamic C++ library, an executable Java JAR file, and a run script.
plots: ${DYNALIGN_SERIAL_FILES} ${PLOT_HANDLER} ${STRUCTURE_IMAGER} ${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend_wrap.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend_wrap.o
	${LINK_PROXY_FOR_JAVA}${DRAWING_PLOTS_LIBRARY} ${DYNALIGN_SERIAL_FILES} ${PLOT_HANDLER} ${STRUCTURE_IMAGER} ${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend_wrap.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend_wrap.o
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_utilities/*.java
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_drawing/source/*/*.java
	cd ${ROOTPATH}; jar cvmf RNAstructure_java_drawing/ManifestPlots.txt exe/RNAstructureDrawingPlots.jar RNAstructure_java_drawing/source/*/*.class RNAstructure_java_utilities/*.class
	sed 's/RNAstructure/RNAstructureDrawingPlots/g' ${ROOTPATH}/exe/RNAstructureScript > ${ROOTPATH}/exe/RNAstructureDrawingPlotsScript
	chmod 755 ${ROOTPATH}/exe/RNAstructureDrawingPlotsScript

# Rule to build the structure drawing program, C++ and Java code.
# This rule builds a dynamic C++ library, an executable Java JAR file, and a run script.
structures: ${HYBRID_FILES} ${STRUCTURE_IMAGER} ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend_wrap.o
	${LINK_PROXY_FOR_JAVA}${DRAWING_STRUCTURES_LIBRARY} ${HYBRID_FILES} ${STRUCTURE_IMAGER} ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.o ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend_wrap.o
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_utilities/*.java
	cd ${ROOTPATH}; ${JXX} ${JXXFLAGS} RNAstructure_java_drawing/source/*/*.java
	cd ${ROOTPATH}; jar cvmf RNAstructure_java_drawing/ManifestStructures.txt exe/RNAstructureDrawingStructures.jar RNAstructure_java_drawing/source/*/*.class RNAstructure_java_utilities/*.class
	sed 's/RNAstructure/RNAstructureDrawingStructures/g' ${ROOTPATH}//exe/RNAstructureScript > ${ROOTPATH}//exe/RNAstructureDrawingStructuresScript
	chmod 755 ${ROOTPATH}//exe/RNAstructureDrawingStructuresScript

# Rule to build SWIG files.
# This rule is almost never necessary for the end user.
swigBuild:
	swig -java -c++ -package RNAstructure_java_drawing.source.proxy -outdir ${ROOTPATH}/RNAstructure_java_drawing/source/proxy ${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.i
	swig -java -c++ -package RNAstructure_java_drawing.source.proxy -outdir ${ROOTPATH}/RNAstructure_java_drawing/source/proxy ${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend.i
	@rm -f ${ROOTPATH}/RNAstructure_java_drawing/source/proxy/*Proxy.java

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

# Rule to clean all possible objects.
clean:
	rm -f ${ROOTPATH}/RNAstructure_java_drawing/source/*/*.class
	rm -f ${ROOTPATH}/RNAstructure_java_utilities/*.class
	rm -f ${PLOT_HANDLER} ${STRUCTURE_IMAGER}
	rm -f ${ROOTPATH}/RNAstructure_java_drawing/SWIG/*.o
	rm -f ${HYBRID_FILES}
	rm -f ${DYNALIGN_SERIAL_FILES}

# Rule to clean all possible objects and executables.
realclean: clean
	rm -f ${ROOTPATH}/exe/${DRAWING_STRUCTURES_LIBRARY} ${ROOTPATH}/exe/RNAstructureDrawingStructures.jar ${ROOTPATH}/exe/RNAstructureDrawingStructuresScript
	rm -f ${ROOTPATH}/exe/${DRAWING_PLOTS_LIBRARY} ${ROOTPATH}/exe/RNAstructureDrawingPlots.jar ${ROOTPATH}/exe/RNAstructureDrawingPlotsScript

