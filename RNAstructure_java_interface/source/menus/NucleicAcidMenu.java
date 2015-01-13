/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.menus;

import RNAstructure_java_interface.source.windows.DynalignFoldWindow;
import RNAstructure_java_interface.source.windows.Efn2Window;
import RNAstructure_java_interface.source.windows.FoldDoubleWindow;
import RNAstructure_java_interface.source.windows.FoldSingleWindow;
import RNAstructure_java_interface.source.windows.FoldSuboptimalWindow;
import RNAstructure_java_interface.source.windows.MaxExpectWindow;
import RNAstructure_java_interface.source.windows.MultilignWindow;
import RNAstructure_java_interface.source.windows.OligoWalkWindow;
import RNAstructure_java_interface.source.windows.PartitionDoubleWindow;
import RNAstructure_java_interface.source.windows.PartitionSingleWindow;
import RNAstructure_java_interface.source.windows.ProbKnotWindow;
import RNAstructure_java_interface.source.windows.PseudoknotWindow;
import RNAstructure_java_interface.source.windows.StochasticWindow;
import RNAstructure_java_interface.source.windows.TurboFoldWindow;

/**
 * A class that creates a menu which holds commands specific for a nucleic
 * acid type.
 * <br><br>
 * Most File menus are identical, but items can be added to or subtracted from
 * them in specific contexts.
 *
 * @author Jessica S. Reuter
 */
public class NucleicAcidMenu
	extends RolloverMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	private NucleicAcidMenu( String acid ) {
		super( acid );
		addMenuItem( "Fold " + acid + " Single Strand",
			"Predict a secondary structure for a single strandusing " + acid +
			" thermodynamics." );
		addMenuItem( "Fold " + acid + " Bimolecular",
			"Predict a secondary structure for two strands using " + acid +
			" thermodynamics." );
		addSeparator();
		addMenuItem( "Partition Function " + acid,
			"Predict base pairing probabilities for all " + acid + " pairs." );
		addMenuItem( "Partition Function " + acid + " Bimolecular",
			"Predict base pairing probabilities for all " + acid + "pairs " +
			"for two strands." );
		addSeparator();
		addMenuItem( "Generate All Suboptimal " + acid + " Structures",
			"Generate all " + acid + " structures within a given increment " +
			"of the lowest free energy structure." );
		addMenuItem( "Stochastic " + acid + " Sampling",
			"Predict " + acid + " structures using stochastic sampling." );
		addMenuItem( "MaxExpect: Predict " + acid + " MEA Structure",
			"Predict the maximum expected accuracy " + acid + " structure." );
		addMenuItem( "ProbKnot: Predict " + acid + " Structures Including " +
			"Pseudoknots",
			"Predict pseudoknots in a " + acid + " structure." );
		addMenuItem( "Efn2 " + acid,
			"Calculate the free energy of a " + acid + " structure." );
		addSeparator();
		addMenuItem( acid + " Dynalign",
			"Find a common secondary structure for two " + acid +
			" sequences." );
		addMenuItem( acid + " Multilign",
			"Predict a secondary structure common to multiple " + acid +
			" sequences." );
		addSeparator();
		addMenuItem( "Break " + acid + " Pseudoknots",
			"Break pseudoknots in a structure, leaving the lowest free " +
			"energy pseudoknot-free structure." );
	}

	/**
	 * Create a DNA menu.
	 *
	 * @return   The DNA menu.
	 */
	public static NucleicAcidMenu createDNA() {
		return new NucleicAcidMenu( "DNA" );
	}

	/**
	 * Create a RNA menu.
	 *
	 * @return   The RNA menu.
	 */
	public static NucleicAcidMenu createRNA() {
		NucleicAcidMenu rna = new NucleicAcidMenu( "RNA" );
		rna.addMenuItem( "RNA OligoWalk", "RNA OligoWalk calculation." );
		rna.insertMenuItem( "RNA TurboFold",
			"Find common structures for multiple sequences using base pair " +
			"probabilities.",
			rna, 15 );
		return rna;
	}

	@Override
	protected void doMenuActions( String command ) {

		// Get the nucleic acid type used for this menu item.
		String type = ( command.contains( "RNA" ) ) ? "RNA" : "DNA";

		// If the command is to break pseudoknots, show a pseudoknot
		// breakage window.
		if( command.startsWith( "Break" ) ) {
			new PseudoknotWindow( type ).viewWindow();
		}

		// If the command is to do a Dynalign calculation, show a Dynalign
		// window.
		else if( command.endsWith( "Dynalign" ) ) {
			new DynalignFoldWindow( type ).viewWindow();
		}

		// If the command is to determine the folding free energy of
		// structures, show an efn2 window.
		else if( command.startsWith( "Efn2" ) ) {
			new Efn2Window( type ).viewWindow();
		}

		// If the command is to do free energy minimization folding, show a
		// folding window.
		else if( command.startsWith( "Fold" ) ) {
			boolean bimol = command.contains( "Bimolecular" );
			if( bimol ) { new FoldDoubleWindow( type ).viewWindow(); }
			else { new FoldSingleWindow( type ).viewWindow(); }
		}

		// If the command is to determine maximum expected accuracy, show a
		// maximum expected accuracy window.
		else if( command.startsWith( "MaxExpect" ) ) {
			new MaxExpectWindow( type ).viewWindow();
		}

		// If the command is to do a Multilign calculation, show a Multilign
		// window.
		else if( command.endsWith( "Multilign" ) ) {
			new MultilignWindow( type ).viewWindow();
		}

		// If the command is to run OligoWalk, show an OligoWalk window.
		else if( command.endsWith( "OligoWalk" ) ) {
			new OligoWalkWindow().viewWindow();
		}

		// If the command is to calculate a partition function, show a
		// partition window.
		else if( command.startsWith( "Partition" ) ) {
			boolean bimol = command.contains( "Bimolecular" );
			if( bimol ) { new PartitionDoubleWindow( type ).viewWindow(); }
			else { new PartitionSingleWindow( type ).viewWindow(); }
		}

		// If the command is to predict pseudoknots, show a pseudoknot
		// prediction window.
		else if( command.startsWith( "ProbKnot" ) ) {
			new ProbKnotWindow( type ).viewWindow();
		}

		// If the command is to perform stochastic sampling, show a stochastic
		// sampling window.
		else if( command.startsWith( "Stochastic" ) ) {
			new StochasticWindow( type ).viewWindow();
		}

		// If the command is to generate all suboptimal structures, show a
		// suboptimal structures window.
		else if( command.contains( "Suboptimal" ) ) {
			new FoldSuboptimalWindow( type ).viewWindow();
		}

		// If the command is to run TurboFold, show a TurboFold window.
		else if( command.contains( "TurboFold" ) ) {
			new TurboFoldWindow().viewWindow();
		}
	}
}
