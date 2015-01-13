/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import RNAstructure_java_interface.source.menus.ConstraintsMenu;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.FileField;
import RNAstructure_java_utilities.FieldPanel.FilePanel;

/**
 * A class responsible for initializing and running the bipartition module,
 * which calculates the partition function for two strands of nucleic acids.
 *
 * @author Jessica S. Reuter
 */
public class PartitionDoubleWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public PartitionDoubleWindow( String acid ) {
		super( acid, acid + " Bimolecular Partition Function" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Save File" button, try to select a
		// partition function save file, and if one was selected set its name.
		if( command.equals( "Save File" ) ) {
			int index = 3;
			String file = FileChooser.doSave( files.getFile( index ) )
				.getPartition();
			if( !file.equals( "" ) ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Sequence File 1" button, try to
		// select a sequence file, and if one was selected set its name.
		else if( command.equals( "Sequence File 1" ) ) {
			String file = FileChooser.doOpen().getSequence();
			if( !file.equals( "" ) ) { files.setFile( 1, file ); }
		}

		// If the action comes from the "Sequence File 2" button, get a second
		// sequence file, initialize a data structure, and create a default
		// output file name.
		else if( command.equals( "Sequence File 2" ) ) {

			// Try to select the second sequence file, and if one was selected
			// set its name. If a sequence file wasn't selected, return.
			String file = FileChooser.doOpen().getSequence();
			if( !file.equals( "" ) ) { files.setFile( 2, file ); }
			else { return; }

			// Get the two input sequence files from the panel.
			String seq1 = files.getFile( 1 );
			String seq2 = files.getFile( 2 );

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildBipartitionDataStructure( seq1, seq2, isRNA );
			if( !moduleInitialized( result ) ) { return; }

			// Build the default output file name and set it in the panel.
			String defaultOut = combineFileNames( seq1, seq2, "pfs" );
			files.setFile( 3, defaultOut );

			// Enable the menus.
			menuBar.enableMenus();
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField seq1 = FileField.createDisabled( "Sequence File 1" );
		FileField seq2 = FileField.createDisabled( "Sequence File 2" );
		FileField pfs = FileField.createEnabled( "Save File" );
		FilePanel files = new FilePanel( this, seq1, seq2, pfs );
		files.setPanelWidth( 400 );
		files.makePanel();

		// Put the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		makeStartButton( 0, 1 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		files.getFile( 2 );
		String pfs = files.getFile( 3 );
		if( files.isError() ) { return; }

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runBipartition( pfs );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// Show the dot plot.
		DrawingWindow drawing = new DrawingWindow( pfs );
		if( drawing.isError() == false ) { drawing.viewWindow(); }
	}

	@Override
	protected RolloverMenu[] setMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
