# Note that the ROOTPATH variable must be defined in a Makefile before this file is included.
# Otherwise, the paths to the various dependencies may not resolve correctly.

##########
## Set variables for the output directory and progress monitoring.
## Use the Java TProgressDialog and ProgressMonitor if libraries are to be integrated with Java.
##########

OBSERVINGBAR = ${ROOTPATH}/src/observingtextprogressbar.o
OUTDIR = ${ROOTPATH}/exe
PROGRESSMONITOR =
TPROGRESSDIR = ${ROOTPATH}/src

ifeq (${JAVA},yes)
	${CXXFLAGS} := -D_JAVA_GUIx
	OBSERVINGBAR =
	OUTDIR = ${ROOTPATH}/RNAstructure_java_interface
	PROGRESSMONITOR = ${ROOTPATH}/RNAstructure_java_interface/SWIG/ProgressMonitor.o
	TPROGRESSDIR = ${ROOTPATH}/RNAstructure_java_interface/SWIG
endif

##########
## Define general convenience macros.
## Note that all the macros in this section are independent, and exist explicitly to make Makefiles clearer.
##########

# The text interface command line parser.
CMD_LINE_PARSER = \
	${ROOTPATH}/src/ParseCommandLine.o

# The text interface configuration file parser.
CONFIG_FILE_PARSER = \
	${ROOTPATH}/src/configfile.o

# The dot plot handler.
PLOT_HANDLER = \
	${ROOTPATH}/src/DotPlotHandler.o

# The utility that creates a structure image comparer.
STRUCTURE_COMPARER = \
	${ROOTPATH}/src/StructureComparedImageHandler.o

# The utility that creates a structure image.
STRUCTURE_IMAGER = \
	${ROOTPATH}/src/StructureImageHandler.o

# The utility that scores a structure against another.
STRUCTURE_SCORER = \
	${ROOTPATH}/src/score.o

##########
## Define file dependency group convenience macros.
##########

# Files unique to PARTS.
PARTS_FILES = \
	${ROOTPATH}/PARTS/src/parts/alignment_priors.o \
	${ROOTPATH}/PARTS/src/parts/array_file_manager.o \
	${ROOTPATH}/PARTS/src/parts/array_mem_manager.o \
	${ROOTPATH}/PARTS/src/parts/map_alignment.o \
	${ROOTPATH}/PARTS/src/parts/map_mhr_info.o \
	${ROOTPATH}/PARTS/src/parts/map_results.o \
	${ROOTPATH}/PARTS/src/parts/map_structures.o \
	${ROOTPATH}/PARTS/src/parts/pf_alignment.o \
	${ROOTPATH}/PARTS/src/parts/pp_results.o \
	${ROOTPATH}/PARTS/src/parts/ppf_cli.o \
	${ROOTPATH}/PARTS/src/parts/ppf_loops.o \
	${ROOTPATH}/PARTS/src/parts/ppf_operators.o \
	${ROOTPATH}/PARTS/src/parts/ppf_progress_bar.o \
	${ROOTPATH}/PARTS/src/parts/ppf_scale.o \
	${ROOTPATH}/PARTS/src/parts/ppf_ss.o \
	${ROOTPATH}/PARTS/src/parts/ppf_tb_stack.o \
	${ROOTPATH}/PARTS/src/parts/ppf_timer.o \
	${ROOTPATH}/PARTS/src/parts/ppf_v_mhe.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w_ext.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w_l.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w_mb.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w_mbl.o \
	${ROOTPATH}/PARTS/src/parts/ppf_w_mhi.o \
	${ROOTPATH}/PARTS/src/parts/process_sequences.o \
	${ROOTPATH}/PARTS/src/parts/single_pf_array.o \
	${ROOTPATH}/PARTS/src/parts/template_pf_array.o \
	${ROOTPATH}/PARTS/src/parts/stoch_tb/stoch_sampled_alignment.o \
	${ROOTPATH}/PARTS/src/parts/stoch_tb/stoch_sampling_math.o \
	${ROOTPATH}/PARTS/src/parts/stoch_tb/stoch_sampled_str_aln_sample_set.o \
	${ROOTPATH}/PARTS/src/parts/stoch_tb/stoch_sampled_structures.o \
	${ROOTPATH}/src/phmm/utils/rng/rng.o \
	${ROOTPATH}/src/phmm/utils/rng/seed_manager.o

# Files for the PHMM calculation group used with both PARTS and TurboFold.
PHMM_FILES = \
	${ROOTPATH}/src/phmm/aln_env_utils.o \
	${ROOTPATH}/src/phmm/p_alignment.o \
	${ROOTPATH}/src/phmm/phmm_aln.o \
	${ROOTPATH}/src/phmm/phmm_array.o \
	${ROOTPATH}/src/phmm/phmm.o \
	${ROOTPATH}/src/phmm/phmm_ml_loops.o \
	${ROOTPATH}/src/phmm/phmm_pp_loops.o \
	${ROOTPATH}/src/phmm/structure/folding_constraints.o \
	${ROOTPATH}/src/phmm/structure/structure_object.o \
	${ROOTPATH}/src/phmm/utils/ansi_string/ansi_string.o \
	${ROOTPATH}/src/phmm/utils/file/utils.o \
	${ROOTPATH}/src/phmm/utils/xmath/linear/linear_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/log/xlog_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/matrix/matrix.o

# Files unique to ShapeKnots.
SHAPEKNOTS_FILES = \
	${ROOTPATH}/ShapeKnots/ShapeKnots_Interface.o \
	${ROOTPATH}/src/ShapeKnots.o \
    ${ROOTPATH}/src/pkHelix.o \
	${ROOTPATH}/src/PseudoParser.o

# TurboFold serial files.
TURBOFOLD_SERIAL_FILES = \
	${ROOTPATH}/TurboFold/TurboFold_object.o \
	${ROOTPATH}/TurboFold/TurboFold_thread.o \
	${RNA_FILES} \
	${PHMM_FILES}

# TurboFold SMP files.
TURBOFOLD_SMP_FILES = \
	${ROOTPATH}/src/phmm/utils/ansi_thread/ansi_thread.o \
	${ROOTPATH}/src/phmm/utils/mutex/ansi_mutex.o \
	${ROOTPATH}/TurboFold/TurboFold_object-smp.o \
	${ROOTPATH}/TurboFold/TurboFold_thread-smp.o \
	${RNA_FILES} \
	${PHMM_FILES}

##########
## Define library dependency group convenience macros.
## Note that this section is not organized alphabetically, but by which macros depend on the ones before them.
##########

# Common files for the RNA library.
RNA_FILES = \
	${ROOTPATH}/RNA_class/RNA.o \
	${ROOTPATH}/RNA_class/thermodynamics.o \
	${ROOTPATH}/src/algorithm.o \
	${ROOTPATH}/src/alltrace.o \
	${ROOTPATH}/src/arrayclass.o \
	${ROOTPATH}/src/dotarray.o \
	${ROOTPATH}/src/draw.o \
	${ROOTPATH}/src/extended_double.o \
	${ROOTPATH}/src/forceclass.o \
	${ROOTPATH}/src/MaxExpect.o \
	${ROOTPATH}/src/MaxExpectStack.o \
	${ROOTPATH}/src/outputconstraints.o \
	${ROOTPATH}/src/pfunction.o \
	${ROOTPATH}/src/probknot.o \
	${ROOTPATH}/src/random.o \
	${ROOTPATH}/src/rna_library.o \
	${ROOTPATH}/src/stackclass.o \
	${ROOTPATH}/src/stackstruct.o \
	${ROOTPATH}/src/stochastic.o \
	${ROOTPATH}/src/structure.o \
	${TPROGRESSDIR}/TProgressDialog.o \
	${PROGRESSMONITOR}


# Common files for the RNA library, for use with OPENMP SMP.
RNA_FILES_SMP = \
	${ROOTPATH}/RNA_class/RNA.o \
	${ROOTPATH}/RNA_class/thermodynamics.o \
	${ROOTPATH}/src/algorithm-smp.o \
	${ROOTPATH}/src/alltrace.o \
	${ROOTPATH}/src/arrayclass.o \
	${ROOTPATH}/src/dotarray.o \
	${ROOTPATH}/src/draw.o \
	${ROOTPATH}/src/extended_double.o \
	${ROOTPATH}/src/forceclass.o \
	${ROOTPATH}/src/MaxExpect.o \
	${ROOTPATH}/src/MaxExpectStack.o \
	${ROOTPATH}/src/outputconstraints.o \
	${ROOTPATH}/src/pfunction-smp.o \
	${ROOTPATH}/src/probknot.o \
	${ROOTPATH}/src/random.o \
	${ROOTPATH}/src/rna_library.o \
	${ROOTPATH}/src/stackclass.o \
	${ROOTPATH}/src/stackstruct.o \
	${ROOTPATH}/src/stochastic-smp.o \
	${ROOTPATH}/src/structure.o \
	${TPROGRESSDIR}/TProgressDialog.o \
	${PROGRESSMONITOR}


RNA_DYNALIGN_II_FILES = \
	${ROOTPATH}/RNA_class/RNA_dynalign_ii.o \
	${ROOTPATH}/RNA_class/thermodynamics.o \
	${ROOTPATH}/src/algorithm_dynalign_ii.o \
	${ROOTPATH}/src/alltrace.o \
	${ROOTPATH}/src/arrayclass.o \
	${ROOTPATH}/src/dotarray.o \
	${ROOTPATH}/src/draw.o \
	${ROOTPATH}/src/extended_double.o \
	${ROOTPATH}/src/forceclass.o \
	${ROOTPATH}/src/MaxExpect.o \
	${ROOTPATH}/src/MaxExpectStack.o \
	${ROOTPATH}/src/outputconstraints.o \
	${ROOTPATH}/src/pfunction.o \
	${ROOTPATH}/src/probknot.o \
	${ROOTPATH}/src/random.o \
	${ROOTPATH}/src/rna_library.o \
	${ROOTPATH}/src/stackclass.o \
	${ROOTPATH}/src/stackstruct.o \
	${ROOTPATH}/src/stochastic.o \
	${ROOTPATH}/src/structure.o \
	${TPROGRESSDIR}/TProgressDialog.o \
	${PROGRESSMONITOR}


# Common files for both versions of the Dynalign library.
DYNALIGN_FILES = \
	${RNA_FILES} \
	${ROOTPATH}/RNA_class/Dynalign_object.o \
	${ROOTPATH}/RNA_class/TwoRNA.o \
	${ROOTPATH}/src/bimol.o \
	${ROOTPATH}/src/dynalignarray.o \
	${ROOTPATH}/src/dynalignheap.o \
	${ROOTPATH}/src/dynalignstackclass.o \
	${ROOTPATH}/src/observable.o \
	${ROOTPATH}/src/observer.o \
	${ROOTPATH}/src/varray.o \
	${ROOTPATH}/src/wendarray.o \
	${ROOTPATH}/src/phmm/aln_env_utils.o \
	${ROOTPATH}/src/phmm/p_alignment.o \
	${ROOTPATH}/src/phmm/phmm.o \
	${ROOTPATH}/src/phmm/phmm_aln.o \
	${ROOTPATH}/src/phmm/phmm_array.o \
	${ROOTPATH}/src/phmm/phmm_ml_loops.o \
	${ROOTPATH}/src/phmm/phmm_pp_loops.o \
	${ROOTPATH}/src/phmm/structure/folding_constraints.o \
	${ROOTPATH}/src/phmm/structure/structure_object.o \
	${ROOTPATH}/src/phmm/utils/ansi_string/ansi_string.o \
	${ROOTPATH}/src/phmm/utils/file/utils.o \
	${ROOTPATH}/src/phmm/utils/xmath/linear/linear_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/matrix/matrix.o \
	${ROOTPATH}/src/phmm/utils/xmath/log/xlog_math.o \
	${CONFIG_FILE_PARSER} \
	${OBSERVINGBAR}


# Common files for both versions of the Dynalign library.
DYNALIGN_II_FILES = \
	${RNA_DYNALIGN_II_FILES} \
	${ROOTPATH}/RNA_class/Dynalign_ii_object.o \
	${ROOTPATH}/RNA_class/TwoRNA.o \
	${ROOTPATH}/src/bimol.o \
	${ROOTPATH}/src/dynalignarray.o \
	${ROOTPATH}/src/dynalignheap.o \
	${ROOTPATH}/src/dynalignstackclass_ii.o \
	${ROOTPATH}/src/observable.o \
	${ROOTPATH}/src/observer.o \
	${ROOTPATH}/src/varray.o \
	${ROOTPATH}/src/wendarray.o \
	${ROOTPATH}/src/phmm/aln_env_utils.o \
	${ROOTPATH}/src/phmm/p_alignment.o \
	${ROOTPATH}/src/phmm/phmm.o \
	${ROOTPATH}/src/phmm/phmm_aln.o \
	${ROOTPATH}/src/phmm/phmm_array.o \
	${ROOTPATH}/src/phmm/phmm_ml_loops.o \
	${ROOTPATH}/src/phmm/phmm_pp_loops.o \
	${ROOTPATH}/src/phmm/structure/folding_constraints.o \
	${ROOTPATH}/src/phmm/structure/structure_object.o \
	${ROOTPATH}/src/phmm/utils/ansi_string/ansi_string.o \
	${ROOTPATH}/src/phmm/utils/file/utils.o \
	${ROOTPATH}/src/phmm/utils/xmath/linear/linear_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/matrix/matrix.o \
	${ROOTPATH}/src/phmm/utils/xmath/log/xlog_math.o \
	${CONFIG_FILE_PARSER} \
	${OBSERVINGBAR}

# Common files for the Dynalign library when compiled serially.
DYNALIGN_SERIAL_FILES = \
	${DYNALIGN_FILES} \
	${ROOTPATH}/src/dynalign.o

# Common files for the Dynalign library when compiling for SMP.
DYNALIGN_SMP_FILES = \
	${DYNALIGN_FILES} \
	${ROOTPATH}/src/dynalign-smp.o \
	${ROOTPATH}/src/rank.o \
	${ROOTPATH}/src/rankconsumer.o \
	${ROOTPATH}/src/rankmanager.o \
	${ROOTPATH}/src/rankproducer.o

# Common files for the Dynalign library when compiled serially.
DYNALIGN_II_SERIAL_FILES = \
	${DYNALIGN_II_FILES} \
	${ROOTPATH}/src/dynalign_ii.o

# Common files for the Dynalign library when compiled serially.
DYNALIGN_II_SMP_FILES = \
	${DYNALIGN_II_FILES} \
	${ROOTPATH}/src/dynalign_ii-smp.o \
	${ROOTPATH}/src/rank.o \
	${ROOTPATH}/src/rankconsumer_ii.o \
	${ROOTPATH}/src/rankmanager.o \
	${ROOTPATH}/src/rankproducer.o

# Common files for the HybridRNA library.
HYBRID_FILES = \
	${RNA_FILES} \
	${ROOTPATH}/RNA_class/HybridRNA.o \
	${ROOTPATH}/RNA_class/TwoRNA.o \
	${ROOTPATH}/src/bimol.o

# Common files for the HybridRNA library, for use with OPENMP SMP.
HYBRID_FILES_SMP = \
	${RNA_FILES_SMP} \
	${ROOTPATH}/RNA_class/HybridRNA.o \
	${ROOTPATH}/RNA_class/TwoRNA.o \
	${ROOTPATH}/src/bimol.o


# Common files for the Oligo library.
OLIGO_FILES = \
	${RNA_FILES} \
	${ROOTPATH}/RNA_class/Oligowalk_object.o \
	${ROOTPATH}/src/alltrace_intermolecular.o \
	${ROOTPATH}/src/intermolecular.o \
	${ROOTPATH}/src/OligoScreenCalc.o \
	${ROOTPATH}/src/pclass.o \
	${ROOTPATH}/src/siRNAfilter.o \
	${ROOTPATH}/src/thermo.o

# Common files for the Oligo library.
OLIGO_FILES_SMP = \
	${RNA_FILES} \
	${ROOTPATH}/RNA_class/Oligowalk_object.o \
	${ROOTPATH}/src/alltrace_intermolecular.o \
	${ROOTPATH}/src/intermolecular.o \
	${ROOTPATH}/src/OligoScreenCalc-smp.o \
	${ROOTPATH}/src/pclass.o \
	${ROOTPATH}/src/siRNAfilter.o \
	${ROOTPATH}/src/thermo.o

# Common files for the phmm library.
PHMM2_FILES = \
	${ROOTPATH}/src/phmm.o

# Common files for the Pseudoknot library.
PSEUDOKNOT_FILES = \
	${ROOTPATH}/src/basepair.o \
	${ROOTPATH}/src/Pseudoknot.o 

# Common files for ProbScan library
PROBSCAN_FILES = \
	${RNA_FILES} \
	${ROOTPATH}/RNA_class/ProbScan.o

# Files for the Java RNAstructure interface (excluding the interface proxies).
# All files are specified explicitly to avoid multiple definitions.
# Updates to this list may be manually necessary if macros change or the Java drawing proxies change.
JAVA_LIBRARY_FILES = \
	${ROOTPATH}/RNA_class/Dynalign_object.o \
	${ROOTPATH}/RNA_class/HybridRNA.o \
	${ROOTPATH}/RNA_class/Multilign_object.o \
	${ROOTPATH}/RNA_class/Oligowalk_object.o \
	${ROOTPATH}/RNA_class/RNA.o \
	${ROOTPATH}/RNA_class/thermodynamics.o \
	${ROOTPATH}/RNA_class/TwoRNA.o \
	${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend.o \
	${ROOTPATH}/RNAstructure_java_drawing/SWIG/DotPlotBackend_wrap.o \
	${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend.o \
	${ROOTPATH}/RNAstructure_java_drawing/SWIG/StructureBackend_wrap.o \
	${ROOTPATH}/src/algorithm.o \
	${ROOTPATH}/src/alltrace.o \
	${ROOTPATH}/src/alltrace_intermolecular.o \
	${ROOTPATH}/src/arrayclass.o \
	${ROOTPATH}/src/bimol.o \
	${ROOTPATH}/src/dotarray.o \
	${ROOTPATH}/src/DotPlotHandler.o \
	${ROOTPATH}/src/draw.o \
	${ROOTPATH}/src/dynalign.o \
	${ROOTPATH}/src/dynalignarray.o \
	${ROOTPATH}/src/dynalignheap.o \
	${ROOTPATH}/src/dynalignstackclass.o \
	${ROOTPATH}/src/extended_double.o \
	${ROOTPATH}/src/forceclass.o \
	${ROOTPATH}/src/intermolecular.o \
	${ROOTPATH}/src/MaxExpect.o \
	${ROOTPATH}/src/MaxExpectStack.o \
	${ROOTPATH}/src/observable.o \
	${ROOTPATH}/src/observer.o \
	${ROOTPATH}/src/OligoScreenCalc.o \
	${ROOTPATH}/src/outputconstraints.o \
	${ROOTPATH}/src/pclass.o \
	${ROOTPATH}/src/pfunction.o \
	${ROOTPATH}/src/probknot.o \
	${ROOTPATH}/src/random.o \
	${ROOTPATH}/src/rna_library.o \
	${ROOTPATH}/src/siRNAfilter.o \
	${ROOTPATH}/src/stackclass.o \
	${ROOTPATH}/src/stackstruct.o \
	${ROOTPATH}/src/stochastic.o \
	${ROOTPATH}/src/structure.o \
	${ROOTPATH}/src/StructureImageHandler.o \
	${ROOTPATH}/src/thermo.o \
	${ROOTPATH}/src/varray.o \
	${ROOTPATH}/src/wendarray.o \
	${ROOTPATH}/src/phmm/aln_env_utils.o \
	${ROOTPATH}/src/phmm/p_alignment.o \
	${ROOTPATH}/src/phmm/phmm_aln.o \
	${ROOTPATH}/src/phmm/phmm_array.o \
	${ROOTPATH}/src/phmm/phmm.o \
	${ROOTPATH}/src/phmm/phmm_ml_loops.o \
	${ROOTPATH}/src/phmm/phmm_pp_loops.o \
	${ROOTPATH}/src/phmm/structure/folding_constraints.o \
	${ROOTPATH}/src/phmm/structure/structure_object.o \
	${ROOTPATH}/src/phmm/utils/ansi_string/ansi_string.o \
	${ROOTPATH}/src/phmm/utils/file/utils.o \
	${ROOTPATH}/src/phmm/utils/xmath/linear/linear_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/log/xlog_math.o \
	${ROOTPATH}/src/phmm/utils/xmath/matrix/matrix.o \
	${ROOTPATH}/TurboFold/TurboFold_object.o \
	${ROOTPATH}/TurboFold/TurboFold_thread.o \
	${TPROGRESSDIR}/TProgressDialog.o \
	${PROGRESSMONITOR}

##########
## Define individual file dependencies.
## Not all files defined in dependency groups above need dependencies here, but most do.
##########

${ROOTPATH}/AllSub/AllSub.o: \
	${ROOTPATH}/AllSub/AllSub.cpp ${ROOTPATH}/AllSub/AllSub.h

${ROOTPATH}/bifold/bifold.o: \
	${ROOTPATH}/bifold/bifold.cpp ${ROOTPATH}/bifold/bifold.h

${ROOTPATH}/CircleCompare/CircleCompare_Interface.o: \
	${ROOTPATH}/CircleCompare/CircleCompare_Interface.cpp ${ROOTPATH}/CircleCompare/CircleCompare_Interface.h \
	${ROOTPATH}/src/StructureComparedImageHandler.cpp ${ROOTPATH}/src/StructureComparedImageHandler.h \
	${ROOTPATH}/src/StructureImageHandler.cpp ${ROOTPATH}/src/StructureImageHandler.h

${ROOTPATH}/dot2ct/dot2ct.o: \
	${ROOTPATH}/dot2ct/dot2ct.cpp ${ROOTPATH}/dot2ct/dot2ct.h

${ROOTPATH}/draw/DrawStructure.o: \
	${ROOTPATH}/draw/DrawStructure.cpp ${ROOTPATH}/draw/DrawStructure.h \
	${ROOTPATH}/src/StructureImageHandler.cpp ${ROOTPATH}/src/StructureImageHandler.h

${ROOTPATH}/DuplexFold/DuplexFold.o: \
	${ROOTPATH}/DuplexFold/DuplexFold.cpp ${ROOTPATH}/DuplexFold/DuplexFold.h

${ROOTPATH}/DynalignDotPlot/DynalignDotPlot.o: \
	${ROOTPATH}/DynalignDotPlot/DynalignDotPlot.cpp ${ROOTPATH}/DynalignDotPlot/DynalignDotPlot.h

${ROOTPATH}/dynalign/dynaligninterface.o: ${ROOTPATH}/dynalign/dynaligninterface.cpp

${ROOTPATH}/dynalign/dynaligninterface-smp.o: ${ROOTPATH}/dynalign/dynaligninterface.cpp
	${COMPILE_SMP} ${ROOTPATH}/dynalign/dynaligninterface.cpp

${ROOTPATH}/dynalign/dynaligninterface_ii.o: ${ROOTPATH}/dynalign/dynaligninterface.cpp
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/dynalign/dynaligninterface.cpp

${ROOTPATH}/dynalign/dynaligninterface_ii-smp.o: ${ROOTPATH}/dynalign/dynaligninterface.cpp
	${COMPILE_DYNALIGN_II_SMP} ${ROOTPATH}/dynalign/dynaligninterface.cpp

${ROOTPATH}/efn2/efn2.o: \
	${ROOTPATH}/efn2/efn2.cpp ${ROOTPATH}/efn2/efn2.h

${ROOTPATH}/efn2/efn2-smp.o: \
	${ROOTPATH}/efn2/efn2.cpp ${ROOTPATH}/efn2/efn2.h
	${CXX} -c ${CXXOPENMPFLAGS} \
	-o ${ROOTPATH}/efn2/efn2-smp.o ${ROOTPATH}/efn2/efn2.cpp
	

${ROOTPATH}/EnergyPlot/EnergyPlot.o: \
	${ROOTPATH}/EnergyPlot/EnergyPlot.cpp ${ROOTPATH}/EnergyPlot/EnergyPlot.h

${ROOTPATH}/EnsembleEnergy/EnsembleEnergy_Interface.o: \
	${ROOTPATH}/EnsembleEnergy/EnsembleEnergy_Interface.cpp ${ROOTPATH}/EnsembleEnergy/EnsembleEnergy_Interface.h

${ROOTPATH}/fold/Fold.o: \
	${ROOTPATH}/fold/Fold.cpp ${ROOTPATH}/fold/Fold.h

${ROOTPATH}/MaxExpect/MaxExpectInterface.o: \
	${ROOTPATH}/MaxExpect/MaxExpectInterface.cpp 	${ROOTPATH}/MaxExpect/MaxExpect.h

${ROOTPATH}/multilign/Multilign_Interface.o: \
	${ROOTPATH}/multilign/Multilign_Interface.cpp ${ROOTPATH}/multilign/Multilign_Interface.h

${ROOTPATH}/multilign/Multilign_Interface-smp.o: \
	${ROOTPATH}/multilign/Multilign_Interface.cpp ${ROOTPATH}/multilign/Multilign_Interface.h
	${COMPILE_SMP} ${ROOTPATH}/multilign/Multilign_Interface.cpp

${ROOTPATH}/Multifind/Multifind_Interface.o: \
	${ROOTPATH}/Multifind/Multifind_Interface.cpp ${ROOTPATH}/Multifind/Multifind_Interface.h
	${COMPILE_SVM} ${ROOTPATH}/Multifind/Multifind_Interface.cpp

${ROOTPATH}/Multifind/Multifind_Interface-smp.o: \
	${ROOTPATH}/Multifind/Multifind_Interface.cpp ${ROOTPATH}/Multifind/Multifind_Interface.h
	${COMPILE_SVM_SMP} ${ROOTPATH}/Multifind/Multifind_Interface.cpp


${ROOTPATH}/napss/napss.o: \
	${ROOTPATH}/napss/napss.cpp \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/configfile.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/oligoscreen/oligoscreen.o: \
	${ROOTPATH}/oligoscreen/oligoscreen.cpp ${ROOTPATH}/oligoscreen/oligoscreen.h

${ROOTPATH}/oligowalk/src/globals.o: \
	${ROOTPATH}/oligowalk/src/globals.cpp ${ROOTPATH}/oligowalk/src/globals.h

${ROOTPATH}/oligowalk/src/oligowalk.o: \
	${ROOTPATH}/oligowalk/src/globals.h \
	${ROOTPATH}/oligowalk/src/oligowalk.cpp \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/alltrace.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/intermolecular.h \
	${ROOTPATH}/src/pclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/siRNAfilter.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/stochastic.h \
	${ROOTPATH}/src/thermo.h \
	${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/pfunction/partition.o: \
	${ROOTPATH}/pfunction/partition.cpp ${ROOTPATH}/pfunction/partition.h

${ROOTPATH}/ProbablePair/ProbablePair.o: \
	${ROOTPATH}/ProbablePair/ProbablePair.cpp ${ROOTPATH}/ProbablePair/ProbablePair.h

${ROOTPATH}/ProbabilityPlot/ProbabilityPlot.o: \
	${ROOTPATH}/ProbabilityPlot/ProbabilityPlot.cpp ${ROOTPATH}/ProbabilityPlot/ProbabilityPlot.h

${ROOTPATH}/ProbKnot/ProbKnot_Interface.o: \
	${ROOTPATH}/ProbKnot/ProbKnot_Interface.cpp ${ROOTPATH}/ProbKnot/ProbKnot_Interface.h

${ROOTPATH}/ProbKnot/ProbScan_Interface.o: \
        ${ROOTPATH}/ProbScan/ProbScan_Interface.cpp ${ROOTPATH}/ProbKnot/ProbScan_Interface.h

${ROOTPATH}/refold/refold.o: \
	${ROOTPATH}/refold/refold.cpp ${ROOTPATH}/refold/refold.h

${ROOTPATH}/RemovePseudoknots/RemovePseudoknots.o: \
	${ROOTPATH}/RemovePseudoknots/RemovePseudoknots.cpp ${ROOTPATH}/RemovePseudoknots/RemovePseudoknots.h

${ROOTPATH}/RNA_class/Dynalign_class.o: RNA_class/Dynalign_class.cpp

${ROOTPATH}/RNA_class/Dynalign_object.o: \
	${ROOTPATH}/RNA_class/Dynalign_object.cpp ${ROOTPATH}/RNA_class/Dynalign_object.h \
	${ROOTPATH}/src/platform.h


${ROOTPATH}/RNA_class/Dynalign_ii_object.o: \
	${ROOTPATH}/RNA_class/Dynalign_object.cpp ${ROOTPATH}/RNA_class/Dynalign_object.h \
	${ROOTPATH}/src/platform.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/RNA_class/Dynalign_object.cpp

${ROOTPATH}/RNA_class/HybridRNA.o: \
	${ROOTPATH}/RNA_class/HybridRNA.cpp ${ROOTPATH}/RNA_class/HybridRNA.h \
	${ROOTPATH}/RNA_class/RNA.cpp ${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.cpp ${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/RNA_class/TwoRNA.cpp ${ROOTPATH}/RNA_class/TwoRNA.h \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/bimol.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h

${ROOTPATH}/RNA_class/HybridRNA_class.o: \
	${ROOTPATH}/RNA_class/HybridRNA_class.cpp \
	${ROOTPATH}/src/structure.h

${ROOTPATH}/RNA_class/Multifind_object.o: \
	${ROOTPATH}/RNA_class/Multifind_object.cpp ${ROOTPATH}/RNA_class/Multifind_object.h
	${COMPILE_SVM_SMP} ${ROOTPATH}/RNA_class/Multifind_object.cpp

${ROOTPATH}/RNA_class/Multilign_object.o: \
	${ROOTPATH}/RNA_class/Dynalign_object.h \
	${ROOTPATH}/RNA_class/Multilign_object.cpp ${ROOTPATH}/RNA_class/Multilign_object.h \
	${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/RNA_class/TwoRNA.h \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/draw.h \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/random.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h


${ROOTPATH}/RNA_class/Multilign_object-Multifind.o: \
	${ROOTPATH}/RNA_class/Dynalign_object.h \
	${ROOTPATH}/RNA_class/Multilign_object.cpp ${ROOTPATH}/RNA_class/Multilign_object.h \
	${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/RNA_class/TwoRNA.h \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/draw.h \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/random.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h
	${COMPILE_MULTIFIND} ${ROOTPATH}/RNA_class/Multilign_object.cpp


${ROOTPATH}/phmm/phmm_interface.o: ${ROOTPATH}/phmm/phmm_interface.cpp \
    ${ROOTPATH}/phmm/phmm_interface.h

${ROOTPATH}/RNA_class/OligoWalk_class.o: ${ROOTPATH}/RNA_class/OligoWalk_class.cpp

${ROOTPATH}/RNA_class/Oligowalk_object.o: \
	${ROOTPATH}/RNA_class/Oligowalk_object.cpp ${ROOTPATH}/RNA_class/Oligowalk_object.h

${ROOTPATH}/RNA_class/RNA.o: \
	${ROOTPATH}/RNA_class/RNA.cpp ${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.cpp ${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/alltrace.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/bimol.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/draw.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/MaxExpect.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/probknot.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/random.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/stochastic.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h

${ROOTPATH}/RNA_class/RNA_dynalign_ii.o: \
	${ROOTPATH}/RNA_class/RNA.cpp ${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.cpp ${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/alltrace.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/bimol.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/draw.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/MaxExpect.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/probknot.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/random.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/stochastic.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/RNA_class/RNA.cpp

${ROOTPATH}/RNA_class/RNA_class.o: \
	${ROOTPATH}/RNA_class/RNA_class.cpp \
	${ROOTPATH}/src/structure.h

${ROOTPATH}/RNA_class/thermodynamics.o: \
	${ROOTPATH}/RNA_class/thermodynamics.cpp ${ROOTPATH}/RNA_class/thermodynamics.h

${ROOTPATH}/RNA_class/TwoRNA.o: \
	${ROOTPATH}/RNA_class/RNA.cpp ${ROOTPATH}/RNA_class/RNA.h \
	${ROOTPATH}/RNA_class/thermodynamics.cpp ${ROOTPATH}/RNA_class/thermodynamics.h \
	${ROOTPATH}/RNA_class/TwoRNA.cpp ${ROOTPATH}/RNA_class/TwoRNA.h \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/bimol.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/pfunction.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h

${ROOTPATH}/RNAstructure_java_interface/SWIG/ProgressMonitor.o: \
	${ROOTPATH}/RNAstructure_java_interface/SWIG/ProgressMonitor.cpp ${ROOTPATH}/RNAstructure_java_interface/SWIG/ProgressMonitor.h

${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator.o: ${ROOTPATH}/RNAstructure_java_interface/SWIG/RNAstructureBackendCalculator.cpp
	${COMPILE_PROXY_FOR_JAVA}

${ROOTPATH}/RNAstructure_java_interface/SWIG/TProgressDialog.o: \
	${ROOTPATH}/RNAstructure_java_interface/SWIG/TProgressDialog.cpp ${ROOTPATH}/RNAstructure_java_interface/SWIG/TProgressDialog.h

${ROOTPATH}/scorer/Scorer_Interface.o: \
	${ROOTPATH}/scorer/Scorer_Interface.cpp ${ROOTPATH}/scorer/Scorer_Interface.h

${ROOTPATH}/src/algorithm.o: \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h

${ROOTPATH}/src/algorithm-smp.o: \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${CXX} -c ${CXXOPENMPFLAGS} \
	-o ${ROOTPATH}/src/algorithm-smp.o ${ROOTPATH}/src/algorithm.cpp 


${ROOTPATH}/src/algorithm_dynalign_ii.o: \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/src/algorithm.cpp

${ROOTPATH}/src/algorithm_instrumented.o: \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${COMPILE_INSTRUMENTED} ${ROOTPATH}/src/algorithm.cpp

${ROOTPATH}/src/alltrace.o: \
	${ROOTPATH}/src/alltrace.cpp ${ROOTPATH}/src/alltrace.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/structure.h

${ROOTPATH}/src/alltrace_intermolecular.o: \
	${ROOTPATH}/src/alltrace_intermolecular.cpp ${ROOTPATH}/src/alltrace_intermolecular.h

${ROOTPATH}/src/arrayclass.o: \
	${ROOTPATH}/src/arrayclass.cpp ${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h

${ROOTPATH}/src/bimol.o: \
	${ROOTPATH}/src/bimol.cpp ${ROOTPATH}/src/bimol.h

${ROOTPATH}/src/configfile.o: \
	${ROOTPATH}/src/configfile.cpp ${ROOTPATH}/src/configfile.h

${ROOTPATH}/src/dotarray.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dotarray.cpp ${ROOTPATH}/src/dotarray.h

${ROOTPATH}/src/DotPlotHandler.o: \
	${ROOTPATH}/src/DotPlotHandler.cpp ${ROOTPATH}/src/DotPlotHandler.h

${ROOTPATH}/src/draw.o: \
	${ROOTPATH}/src/draw.cpp ${ROOTPATH}/src/draw.h \
	${ROOTPATH}/src/substructure.cpp

${ROOTPATH}/src/dynalign.o: \
    ${ROOTPATH}/src/algorithm.cpp

${ROOTPATH}/src/dynalign.o: \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.cpp ${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/dynalignstackclass.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h \
	${TPROGRESSDIR}/TProgressDialog.h

${ROOTPATH}/src/dynalign-smp.o: \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.cpp ${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/dynalignstackclass.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/observingtextprogressbar.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rankconsumer.h \
	${ROOTPATH}/src/rankmanager.h \
	${ROOTPATH}/src/rankproducer.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h
	${COMPILE_SMP} ${ROOTPATH}/src/dynalign.cpp

${ROOTPATH}/src/dynalign_ii.o: \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.cpp ${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/dynalignstackclass.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/src/dynalign.cpp

${ROOTPATH}/src/dynalign_ii-smp.o: \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.cpp ${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/dynalignheap.h \
	${ROOTPATH}/src/dynalignstackclass.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h \
	${TPROGRESSDIR}/TProgressDialog.h
	${COMPILE_DYNALIGN_II_SMP} ${ROOTPATH}/src/dynalign.cpp

${ROOTPATH}/src/dynalignarray.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.cpp ${ROOTPATH}/src/dynalignarray.h

${ROOTPATH}/src/dynalignheap.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalignheap.cpp ${ROOTPATH}/src/dynalignheap.h

${ROOTPATH}/src/dynalignstackclass.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalignstackclass.cpp ${ROOTPATH}/src/dynalignstackclass.h


${ROOTPATH}/src/dynalignstackclass_ii.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalignstackclass.cpp ${ROOTPATH}/src/dynalignstackclass.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/src/dynalignstackclass.cpp

${ROOTPATH}/src/extended_double.o: \
	${ROOTPATH}/src/extended_double.cpp ${ROOTPATH}/src/extended_double.h

${ROOTPATH}/src/forceclass.o: \
	${ROOTPATH}/src/forceclass.cpp ${ROOTPATH}/src/forceclass.h

${ROOTPATH}/src/intermolecular.o: \
	${ROOTPATH}/src/intermolecular.cpp ${ROOTPATH}/src/intermolecular.h \
	${ROOTPATH}/src/siRNAfilter.cpp ${ROOTPATH}/src/siRNAfilter.h

${ROOTPATH}/src/MaxExpect.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/MaxExpect.cpp ${ROOTPATH}/src/MaxExpect.h

${ROOTPATH}/src/MaxExpectStack.o: \
	${ROOTPATH}/src/MaxExpectStack.cpp ${ROOTPATH}/src/MaxExpectStack.h

${ROOTPATH}/src/observable.o: \
	${ROOTPATH}/src/observable.cpp ${ROOTPATH}/src/observable.h \
	${ROOTPATH}/src/observer.h

${ROOTPATH}/src/observer.o: \
	${ROOTPATH}/src/observer.cpp ${ROOTPATH}/src/observer.h

${ROOTPATH}/src/observingtextprogressbar.o: \
	${ROOTPATH}/src/observingtextprogressbar.cpp ${ROOTPATH}/src/observingtextprogressbar.h \
	${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/src/OligoScreenCalc.o: \
	${ROOTPATH}/src/OligoScreenCalc.cpp ${ROOTPATH}/src/OligoScreenCalc.h

${ROOTPATH}/src/OligoScreenCalc-smp.o: \
	${ROOTPATH}/src/OligoScreenCalc.cpp ${ROOTPATH}/src/OligoScreenCalc.h
	${CXX} -c ${CXXOPENMPFLAGS} \
	-o ${ROOTPATH}/src/OligoScreenCalc-smp.o ${ROOTPATH}/src/OligoScreenCalc.cpp 


${ROOTPATH}/src/outputconstraints.o: \
	${ROOTPATH}/src/outputconstraints.cpp ${ROOTPATH}/src/outputconstraints.h

${ROOTPATH}/src/ParseCommandLine.o: \
	${ROOTPATH}/src/ParseCommandLine.cpp ${ROOTPATH}/src/ParseCommandLine.h

${ROOTPATH}/src/pclass.o: \
	${ROOTPATH}/src/pclass.cpp ${ROOTPATH}/src/pclass.h

${ROOTPATH}/src/probknot.o: \
	${ROOTPATH}/src/probknot.cpp ${ROOTPATH}/src/probknot.h

${ROOTPATH}/src/pfunction.o: \
	${ROOTPATH}/src/pfunction.cpp ${ROOTPATH}/src/pfunction.h ${ROOTPATH}/src/boltzmann.h \
	${ROOTPATH}/src/algorithm.h ${ROOTPATH}/src/structure.h  

${ROOTPATH}/src/pfunction-smp.o: \
	${ROOTPATH}/src/pfunction.cpp ${ROOTPATH}/src/pfunction.h ${ROOTPATH}/src/boltzmann.h \
	${ROOTPATH}/src/algorithm.h ${ROOTPATH}/src/structure.h  
	${CXX} -c ${CXXOPENMPFLAGS} \
	-o ${ROOTPATH}/src/pfunction-smp.o ${ROOTPATH}/src/pfunction.cpp

${ROOTPATH}/src/phmm.o: \
    ${ROOTPATH}/src/phmm.cpp ${ROOTPATH}/src/phmm.h

${ROOTPATH}/src/random.o: \
	${ROOTPATH}/src/random.cpp ${ROOTPATH}/src/random.h

${ROOTPATH}/src/rank.o: \
	${ROOTPATH}/src/rank.cpp ${ROOTPATH}/src/rank.h \
	${ROOTPATH}/src/workslice.h \
	${ROOTPATH}/src/workunit.h

${ROOTPATH}/src/rankconsumer.o: \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/rankconsumer.cpp ${ROOTPATH}/src/rankconsumer.h \
	${ROOTPATH}/src/rankmanager.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h \
	${ROOTPATH}/src/workslice.h

${ROOTPATH}/src/rankconsumer_ii.o: \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/dynalignarray.h \
	${ROOTPATH}/src/rankconsumer.cpp ${ROOTPATH}/src/rankconsumer.h \
	${ROOTPATH}/src/rankmanager.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/varray.h \
	${ROOTPATH}/src/wendarray.h \
	${ROOTPATH}/src/workslice.h
	${COMPILE_DYNALIGN_II} ${ROOTPATH}/src/rankconsumer.cpp

${ROOTPATH}/src/rankmanager.o: \
	${ROOTPATH}/src/observable.h \
	${ROOTPATH}/src/rank.h \
	${ROOTPATH}/src/rankmanager.cpp ${ROOTPATH}/src/rankmanager.h \
	${ROOTPATH}/src/TProgressDialog.h \
	${ROOTPATH}/src/workslice.h

${ROOTPATH}/src/rankproducer.o: \
	${ROOTPATH}/src/rank.h \
	${ROOTPATH}/src/rankmanager.h \
	${ROOTPATH}/src/rankproducer.cpp ${ROOTPATH}/src/rankproducer.h \
	${ROOTPATH}/src/workunit.h

${ROOTPATH}/src/rna_library.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.cpp ${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h

${ROOTPATH}/src/pkHelix.o: \
	${ROOTPATH}/src/pkHelix.cpp ${ROOTPATH}/src/pkHelix.h

${ROOTPATH}/src/basepair.o: \
	${ROOTPATH}/src/basepair.cpp ${ROOTPATH}/src/basepair.h

${ROOTPATH}/src/Pseudoknot.o: \
	${ROOTPATH}/src/Pseudoknot.cpp ${ROOTPATH}/src/Pseudoknot.h

${ROOTPATH}/src/PseudoParser.o: \
	${ROOTPATH}/src/PseudoParser.cpp ${ROOTPATH}/src/PseudoParser.h

${ROOTPATH}/src/score.o: \
	${ROOTPATH}/src/score.cpp ${ROOTPATH}/src/score.h

${ROOTPATH}/src/siRNAfilter.o: \
	${ROOTPATH}/src/siRNAfilter.cpp ${ROOTPATH}/src/siRNAfilter.h

${ROOTPATH}/src/stackclass.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/stackclass.cpp ${ROOTPATH}/src/stackclass.h

${ROOTPATH}/src/stackstruct.o: \
	${ROOTPATH}/src/stackstruct.cpp ${ROOTPATH}/src/stackstruct.h

${ROOTPATH}/src/stochastic.o: \
	${ROOTPATH}/src/stochastic.cpp ${ROOTPATH}/src/stochastic.h

${ROOTPATH}/src/stochastic-smp.o: \
	${ROOTPATH}/src/stochastic.cpp ${ROOTPATH}/src/stochastic.h
	${CXX} -c ${CXXOPENMPFLAGS} \
	-o ${ROOTPATH}/src/stochastic-smp.o ${ROOTPATH}/src/stochastic.cpp 

${ROOTPATH}/src/structure.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/structure.cpp ${ROOTPATH}/src/structure.h

${ROOTPATH}/src/StructureComparedImageHandler.o: \
	${ROOTPATH}/src/StructureComparedImageHandler.cpp ${ROOTPATH}/src/StructureComparedImageHandler.h \
	${ROOTPATH}/src/StructureImageHandler.cpp ${ROOTPATH}/src/StructureImageHandler.h

${ROOTPATH}/src/StructureImageHandler.o: \
	${ROOTPATH}/src/StructureImageHandler.cpp ${ROOTPATH}/src/StructureImageHandler.h

${ROOTPATH}/src/thermo.o: \
	${ROOTPATH}/src/thermo.cpp ${ROOTPATH}/src/thermo.h

${ROOTPATH}/src/TProgressDialog.o: \
	${ROOTPATH}/src/TProgressDialog.cpp ${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/src/varray.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/varray.cpp ${ROOTPATH}/src/varray.h

${ROOTPATH}/src/wendarray.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/wendarray.cpp ${ROOTPATH}/src/wendarray.h

${ROOTPATH}/stochastic/stochastic.o: \
	${ROOTPATH}/stochastic/stochastic.cpp ${ROOTPATH}/stochastic/stochastic.h

${ROOTPATH}/ShapeKnots/ShapeKnots_Interface.o: \
    ${ROOTPATH}/ShapeKnots/ShapeKnots_Interface.cpp ${ROOTPATH}/ShapeKnots/ShapeKnots_Interface.h \
    ${ROOTPATH}/src/ShapeKnots.h 

${ROOTPATH}/src/ShapeKnots.o: \
    ${ROOTPATH}/src/ShapeKnots.cpp ${ROOTPATH}/src/ShapeKnots.h \
	${ROOTPATH}/src/pkHelix.h \
	${ROOTPATH}/src/PseudoParser.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/ParseCommandLine.h \
	${ROOTPATH}/RNA_class/RNA.h 

${ROOTPATH}/TurboFold/TurboFold_Interface-smp.o:
	${COMPILE_SMP} ${ROOTPATH}/TurboFold/TurboFold_Interface.cpp

${ROOTPATH}/TurboFold/TurboFold_object-smp.o:
	${COMPILE_SMP} ${ROOTPATH}/TurboFold/TurboFold_object.cpp

${ROOTPATH}/TurboFold/TurboFold_thread-smp.o:
	${COMPILE_SMP} ${ROOTPATH}/TurboFold/TurboFold_thread.cpp

