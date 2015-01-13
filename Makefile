##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'Use "make all" to create all executables.'
	@echo 'Use "make serial" to create all the serial executables.'
	@echo 'Use "make SMP" to create all available SMP parallel executables.'
	@echo 'Use "make AllSub" to create executable "AllSub."'
	@echo 'Use "make bifold" to create executable "bifold."'
	@echo 'Use "make bipartition" to create executable "bipartition."'
	@echo 'Use "make CircleCompare" to create executable "CircleCompare."'
	@echo 'Use "make ct2dot" to create executable "ct2dot."'
	@echo 'Use "make dot2ct" to create executable "dot2ct."'
	@echo 'Use "make draw" to create executable "draw."'
	@echo 'Use "make dynalign" to create exectuable "dynalign."'
	@echo 'Use "make dynalign-smp" to create executable "dynalign-smp."'
	@echo 'Use "make dynalign_ii" to create exectuable "dynalign_ii."'
	@echo 'Use "make dynalign_ii-smp" to create executable "dynalign_ii-smp."'
	@echo 'Use "make DuplexFold" to create executable "DuplexFold."'
	@echo 'Use "make DynalignDotPlot" to create executable "DynalignDotPlot."'
	@echo 'Use "make efn2" to create executable "efn2."'
	@echo 'Use "make EnergyPlot" to create executable "EnergyPlot."'
	@echo 'Use "make EnsembleEnergy" to create executable "EnsembleEnergy."'
	@echo 'Use "make Fold" to create executable "Fold."'
	@echo 'Use "make fold-cuda" to create executable "fold-cuda."'
	@echo 'Use "make MaxExpect" to create executable "MaxExpect."'
	@echo 'Use "make multilign" to create executable "multilign."'
	@echo 'Use "make multilign-smp" to create executable "multilign-smp."'
	@echo 'Use "make NAPSS" to create the executable "NAPSS."'
	@echo 'Use "make oligoscreen" to create executable "oligoscreen."'
	@echo 'Use "make OligoWalk" to create executable "OligoWalk."'
	@echo 'Use "make partition" to create executable "partition."'
	@echo 'Use "make partition-cuda" to create executable "partition-cuda."'
	@echo 'Use "make PARTS" to create executable "PARTS."'
	@echo 'Use "make phmm" to create executable "phmm."'	
	@echo 'Use "make ProbabilityPlot" to create executable "ProbabilityPlot."'
	@echo 'Use "make ProbablePair" to create executable "ProbablePair."'
	@echo 'Use "make ProbKnot" to create executable "ProbKnot."'
	@echo 'Use "make ProbScan" to create executable "ProbScan."'
	@echo 'Use "make refold" to create executable "refold."'
	@echo 'Use "make RemovePseudoknots" to create executable "RemovePseudoknots."'
	@echo 'Use "make scorer" to create executable "scorer."'
	@echo 'Use "make ShapeKnots" to create executable "ShapeKnots."'	
	@echo 'Use "make stochastic" to create executable "stochastic."'
	@echo 'Use "make TurboFold" to create executable "TurboFold."'
	@echo 'Use "make TurboFold-smp" to create executable "TurboFold-smp."'

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=.
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

##########
## Define repository build actions.
##########

# Make all the executables at once.
all:
	make serial;
	make SMP;

# Make all the serial executables.
serial:
	@echo "Building of all RNAstructure serial programs started."
	@echo
	make AccessFold;
	make AllSub;
	make bifold;
	make bipartition;
	make CircleCompare;
	make ct2dot;
	make dot2ct;
	make draw;
	make DuplexFold;
	make dynalign;
	make dynalign_ii;
	make DynalignDotPlot;
	make efn2;
	make EnergyPlot;
	make EnsembleEnergy;
	make Fold;
	make MaxExpect;
	make multilign;
	make NAPSS;
	make oligoscreen;
	make OligoWalk;
	make partition;
	make PARTS;
	make phmm;
	make ProbabilityPlot;
	make ProbablePair;
	make ProbKnot;
	make refold;
	make RemovePseudoknots;
	make scorer;
	make ShapeKnots;
	make stochastic;
	make TurboFold;
	@echo
	@echo "Building of the serial RNAstructure programs finished."

# Make all the parallel executables at once.
SMP:
	@echo "Building of all RNAstructure SMP programs started."
	@echo
	make bifold-smp;
	make bipartition-smp;
	make dynalign-smp;
	make dynalign_ii-smp;
	make efn2-smp
	make Fold-smp;
	make multilign-smp;
	make partition-smp;
	make oligoscreen-smp;
	make stochastic-smp;
	make TurboFold-smp;
	@echo
	@echo "Building of the parallel RNAstructure programs finished."

# Copy the executables to the /usr/local directory.
install:
	cp -r exe/ /usr/local/RNAstructure

##########
## Define targets.
##########

# Build the AccessFold text interface.
AccessFold: exe/AccessFold
exe/AccessFold: AccessFold/AccessFold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}
	${LINK} AccessFold/AccessFold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}


# Build the AllSub text interface.
AllSub: exe/AllSub
exe/AllSub: AllSub/AllSub.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} AllSub/AllSub.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the bifold text interface.
bifold: exe/bifold
exe/bifold: bifold/bifold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}
	${LINK} bifold/bifold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}

# Build the bifold-smp text interface.
bifold-smp: exe/bifold-smp
exe/bifold-smp: bifold/bifold.o ${CMD_LINE_PARSER} ${HYBRID_FILES_SMP}
	${LINKSMP} bifold/bifold.o ${CMD_LINE_PARSER} ${HYBRID_FILES_SMP}

# Build the bipartition text interface.
bipartition: exe/bipartition
exe/bipartition: bipartition/bipartition.o ${CMD_LINE_PARSER} ${HYBRID_FILES}
	${LINK} bipartition/bipartition.o ${CMD_LINE_PARSER} ${HYBRID_FILES}

# Build the bipartition-smp text interface.
bipartition-smp: exe/bipartition-smp
exe/bipartition-smp: bipartition/bipartition.o ${CMD_LINE_PARSER} ${HYBRID_FILES_SMP}
	${LINKSMP} bipartition/bipartition.o ${CMD_LINE_PARSER} ${HYBRID_FILES_SMP}

# Build the CircleCompare text interface.
CircleCompare: exe/CircleCompare
exe/CircleCompare: CircleCompare/CircleCompare_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${STRUCTURE_IMAGER} ${STRUCTURE_COMPARER} ${RNA_FILES}
	${LINK} CircleCompare/CircleCompare_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${STRUCTURE_IMAGER} ${STRUCTURE_COMPARER} ${RNA_FILES}

# Build the ct2dot text interface.
ct2dot: exe/ct2dot
exe/ct2dot: ct2dot/ct2dot.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} ct2dot/ct2dot.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the dot2ct text interface.
dot2ct: exe/dot2ct
exe/dot2ct: dot2ct/dot2ct.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} dot2ct/dot2ct.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the draw text interface.
draw: exe/draw
exe/draw: draw/DrawStructure.o ${CMD_LINE_PARSER} ${STRUCTURE_IMAGER} ${RNA_FILES}
	${LINK} draw/DrawStructure.o ${CMD_LINE_PARSER} ${STRUCTURE_IMAGER} ${RNA_FILES}

# Build the DuplexFold text interface.
DuplexFold: exe/DuplexFold
exe/DuplexFold: DuplexFold/DuplexFold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}
	${LINK} DuplexFold/DuplexFold.o ${CMD_LINE_PARSER} ${HYBRID_FILES}

# Build the dynalign serial text interface.
dynalign: exe/dynalign
exe/dynalign: dynalign/dynaligninterface.o ${DYNALIGN_SERIAL_FILES}
	${LINK} dynalign/dynaligninterface.o ${DYNALIGN_SERIAL_FILES}

# Build the dynalign SMP text interface.
dynalign-smp: exe/dynalign-smp
exe/dynalign-smp: dynalign/dynaligninterface-smp.o ${DYNALIGN_SMP_FILES}
	${LINK} dynalign/dynaligninterface-smp.o ${DYNALIGN_SMP_FILES} -lpthread

# Build the dynalign serial text interface.
dynalign_ii: exe/dynalign_ii
exe/dynalign_ii: dynalign/dynaligninterface_ii.o ${DYNALIGN_II_SERIAL_FILES}
	${LINK} dynalign/dynaligninterface_ii.o ${DYNALIGN_II_SERIAL_FILES}

dynalign_ii-smp: exe/dynalign_ii-smp
exe/dynalign_ii-smp: dynalign/dynaligninterface_ii-smp.o ${DYNALIGN_II_SMP_FILES}
	${LINK} dynalign/dynaligninterface_ii-smp.o ${DYNALIGN_II_SMP_FILES} -lpthread

# Build the DynalignDotPlot text interface.
DynalignDotPlot: exe/DynalignDotPlot
exe/DynalignDotPlot: DynalignDotPlot/DynalignDotPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${DYNALIGN_SERIAL_FILES}
	${LINK} DynalignDotPlot/DynalignDotPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${DYNALIGN_SERIAL_FILES}

# Build the efn2 text interface.
efn2: exe/efn2
exe/efn2: efn2/efn2.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} efn2/efn2.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the efn2-smp text interface.
efn2-smp: exe/efn2-smp
exe/efn2-smp: efn2/efn2-smp.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINKSMP} efn2/efn2-smp.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the EnergyPlot text interface.
EnergyPlot: exe/EnergyPlot
exe/EnergyPlot: EnergyPlot/EnergyPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${RNA_FILES}
	${LINK} EnergyPlot/EnergyPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${RNA_FILES}

# Build the EnsembleEnergy text interface.
EnsembleEnergy: exe/EnsembleEnergy
exe/EnsembleEnergy: EnsembleEnergy/EnsembleEnergy_Interface.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} EnsembleEnergy/EnsembleEnergy_Interface.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the Fold text interface.
Fold: exe/Fold
exe/Fold: fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the Fold-smp text interface.
Fold-smp: exe/Fold-smp
exe/Fold-smp: fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES_SMP}
	${LINKSMP} fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES_SMP}

# Build the MaxExpect text interface.
MaxExpect: exe/MaxExpect
exe/MaxExpect: MaxExpect/MaxExpectInterface.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} MaxExpect/MaxExpectInterface.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the multilign serial text interface.
multilign: exe/multilign
exe/multilign: multilign/Multilign_Interface.o RNA_class/Multilign_object.o ${CMD_LINE_PARSER} ${DYNALIGN_SERIAL_FILES}
	${LINK} multilign/Multilign_Interface.o RNA_class/Multilign_object.o ${CMD_LINE_PARSER} ${DYNALIGN_SERIAL_FILES}

# Build the multilign SMP text interface.
multilign-smp: exe/multilign-smp
exe/multilign-smp: multilign/Multilign_Interface.o RNA_class/Multilign_object.o ${CMD_LINE_PARSER} ${DYNALIGN_SMP_FILES}
	${LINK} multilign/Multilign_Interface.o RNA_class/Multilign_object.o ${CMD_LINE_PARSER} ${DYNALIGN_SMP_FILES} -lpthread

Multifind: Multifind/Multifind_Interface.o RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o ${CMD_LINE_PARSER} ${DYNALIGN_SERIAL_FILES} ${INCLUDESVM}/svm.o 
#	cd ${INCLUDESVM};make svm.o;
	${LINK_IN_EXE_DIR} Multifind/Multifind_Interface.o  RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o ${CMD_LINE_PARSER} ${DYNALIGN_SERIAL_FILES} ${INCLUDESVM}/svm.o 

Multifind-smp: Multifind/Multifind_Interface-smp.o RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o ${CMD_LINE_PARSER} ${DYNALIGN_SMP_FILES} ${INCLUDESVM}/svm.o 
#	cd ${INCLUDESVM};make svm.o;
	${LINK_IN_EXE_DIR} Multifind/Multifind_Interface-smp.o  RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o ${CMD_LINE_PARSER} ${DYNALIGN_SMP_FILES} ${INCLUDESVM}/svm.o -lpthread

# Build the NAPSS text interface.
NAPSS: exe/NAPSS
exe/NAPSS: napss/napss.o src/algorithm_instrumented.o ${CMD_LINE_PARSER} ${CONFIG_FILE_PARSER} ${RNA_FILES}
	${LINK} napss/napss.o src/algorithm_instrumented.o ${CMD_LINE_PARSER} ${CONFIG_FILE_PARSER} ${RNA_FILES}

# Build the oligoscreen text interface.
oligoscreen: exe/oligoscreen
exe/oligoscreen: oligoscreen/oligoscreen.o ${CMD_LINE_PARSER} ${OLIGO_FILES}
	${LINK} oligoscreen/oligoscreen.o ${CMD_LINE_PARSER} ${OLIGO_FILES}

# Build the oligoscreen text interface.
oligoscreen-smp: exe/oligoscreen-smp
exe/oligoscreen-smp: oligoscreen/oligoscreen.o ${CMD_LINE_PARSER} ${OLIGO_FILES_SMP}
	${LINKSMP} oligoscreen/oligoscreen.o ${CMD_LINE_PARSER} ${OLIGO_FILES_SMP}


# Build the OligoWalk text interface.
OligoWalk: exe/OligoWalk
exe/OligoWalk: oligowalk/src/oligowalk.o oligowalk/src/globals.o ${OLIGO_FILES}
	${LINK} oligowalk/src/oligowalk.o oligowalk/src/globals.o ${OLIGO_FILES}

# Build the partition text interface.
partition: exe/partition
exe/partition: pfunction/partition.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} pfunction/partition.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build partition-cuda
partition-cuda:
	cd partition-smp ; make ../exe/partition-cuda

# Build fold-cuda
fold-cuda:
	cd fold-smp ; make ../exe/fold-cuda

# Build the partition-smp text interface for SMP.
partition-smp: exe/partition-smp
exe/partition-smp: pfunction/partition.o ${CMD_LINE_PARSER} ${RNA_FILES_SMP}
	${LINKSMP} pfunction/partition.o ${CMD_LINE_PARSER} ${RNA_FILES_SMP}

# Build the PARTS text interface.
# This target is unique in RNAstructure and is used to avoid name conflicts.
PARTS: PARTS-recursive
PARTS-recursive:
	cd PARTS; make PARTS

# Build the phmm text interface.
phmm: exe/phmm
exe/phmm: phmm/phmm_interface.o ${DYNALIGN_SERIAL_FILES} ${PHMM2} ${CMD_LINE_PARSER}
	${LINK} phmm/phmm_interface.o ${DYNALIGN_SERIAL_FILES} ${PHMM2} ${CMD_LINE_PARSER}

# Build the ProbabilityPlot text interface.
ProbabilityPlot: exe/ProbabilityPlot
exe/ProbabilityPlot: ProbabilityPlot/ProbabilityPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${RNA_FILES}
	${LINK} ProbabilityPlot/ProbabilityPlot.o ${CMD_LINE_PARSER} ${PLOT_HANDLER} ${RNA_FILES}

# Build the ProbablePair text interface.
ProbablePair: exe/ProbablePair
exe/ProbablePair: ProbablePair/ProbablePair.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} ProbablePair/ProbablePair.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the ProbKnot text interface.
ProbKnot: exe/ProbKnot
exe/ProbKnot: ProbKnot/ProbKnot_Interface.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} ProbKnot/ProbKnot_Interface.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the ProbScan text interface.
ProbScan: exe/ProbScan
exe/ProbScan: ProbScan/ProbScan_Interface.o ${CMD_LINE_PARSER} ${PROBSCAN_FILES}
	${LINK} ProbScan/ProbScan_Interface.o ${CMD_LINE_PARSER} ${PROBSCAN_FILES}

# Build the refold text interface.
refold: exe/refold
exe/refold: refold/refold.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} refold/refold.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the RemovePseudoknots text interface.
RemovePseudoknots: exe/RemovePseudoknots
exe/RemovePseudoknots: RemovePseudoknots/RemovePseudoknots.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} RemovePseudoknots/RemovePseudoknots.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the scorer interface.
scorer: exe/scorer
exe/scorer: scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	${LINK} scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}

# Build the ShapeKnots text interface.
ShapeKnots: exe/ShapeKnots
exe/ShapeKnots: ${CMD_LINE_PARSER} ${RNA_FILES} ${PSEUDOKNOT_FILES} ${SHAPEKNOTS_FILES}
	${LINK} ${CMD_LINE_PARSER} ${RNA_FILES} ${PSEUDOKNOT_FILES} ${SHAPEKNOTS_FILES}


# Build the stochastic text interface.
stochastic: exe/stochastic
exe/stochastic: stochastic/stochastic.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} stochastic/stochastic.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the stochastic text interface.
stochastic-smp: exe/stochastic-smp
exe/stochastic-smp: stochastic/stochastic.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINKSMP} stochastic/stochastic.o ${CMD_LINE_PARSER} ${RNA_FILES_SMP}


# Build the TurboFold serial text interface.
# This target is unique in RNAstructure and is used to avoid name conflicts.
TurboFold: TurboFold-recursive
TurboFold-recursive:
	cd TurboFold; make TurboFold

# Build the TurboFold SMP text interface.
TurboFold-smp:
	cd TurboFold; make TurboFold-smp

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

# Remove object files and any temporary files from building.
clean:
	find . -depth -name '*~' -delete
	find . -depth -name '*.o' -delete
	find . -depth -name '*.class' -delete

# Remove object files and executables.
realclean: clean
	rm -f RNA_class/*_class
	find exe -maxdepth 1 -type f ! -name RNAstructureScript ! -name RNAstructure.bat -delete
	find exe -type d ! -name exe ! -name CVS | xargs rm -rf

