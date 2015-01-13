##########
## Set Library Names
## Options for Windows, Mac, and Linux (default)
##########

ifeq (${OPSYSTEM},Windows)
	DRAWING_LIBRARY = Drawing.dll
	DRAWING_PLOTS_LIBRARY = DrawingPlots.dll
	DRAWING_STRUCTURES_LIBRARY = DrawingStructures.dll
	DYNALIGN_LIBRARY = Dynalign.dll
	DYNALIGN_SMP_LIBRARY = Dynalign_SMP.dll
	HYBRID_RNA_LIBRARY = HybridRNA.dll
	MULTILIGN_LIBRARY = Multilign.dll
	MULTILIGN_SMP_LIBRARY = Multilign_SMP.dll
	OLIGO_LIBRARY = Oligo.dll
	PARTS_LIBRARY = PARTS.dll
	RNA_LIBRARY = RNA.dll
	RNASTRUCTURE_LIBRARY = RNAstructure_GUI.dll
	TURBOFOLD_LIBRARY = TurboFold.dll
	TURBOFOLD_SMP_LIBRARY = TurboFold_SMP.dll
else ifeq (${OPSYSTEM},Mac)
	DRAWING_LIBRARY = libDrawing.dylib
	DRAWING_PLOTS_LIBRARY = libDrawingPlots.dylib
	DRAWING_STRUCTURES_LIBRARY = libDrawingStructures.dylib
	DYNALIGN_LIBRARY = libDynalign.dylib
	DYNALIGN_SMP_LIBRARY = libDynalign_SMP.dylib
	HYBRID_RNA_LIBRARY = libHybridRNA.dylib
	MULTILIGN_LIBRARY = libMultilign.dylib
	MULTILIGN_SMP_LIBRARY = libMultilign_SMP.dylib
	OLIGO_LIBRARY = libOligo.dylib
	PARTS_LIBRARY = libPARTS.dylib
	RNA_LIBRARY = libRNA.dylib
	RNASTRUCTURE_LIBRARY = libRNAstructure_GUI.dylib
	TURBOFOLD_LIBRARY = libTurboFold.dylib
	TURBOFOLD_SMP_LIBRARY = libTurboFold_SMP.dylib
else
	DRAWING_LIBRARY = libDrawing.so
	DRAWING_PLOTS_LIBRARY = libDrawingPlots.so
	DRAWING_STRUCTURES_LIBRARY = libDrawingStructures.so
	DYNALIGN_LIBRARY = libDynalign.so
	DYNALIGN_SMP_LIBRARY = libDynalign_SMP.so
	HYBRID_RNA_LIBRARY = libHybridRNA.so
	MULTILIGN_LIBRARY = libMultilign.so
	MULTILIGN_SMP_LIBRARY = libMultilign_SMP.so
	OLIGO_LIBRARY = libOligo.so
	PARTS_LIBRARY = libPARTS.so
	RNA_LIBRARY = libRNA.so
	RNASTRUCTURE_LIBRARY = libRNAstructure_GUI.so
	TURBOFOLD_LIBRARY = libTurboFold.so
	TURBOFOLD_SMP_LIBRARY = libTurboFold_SMP.so
endif

