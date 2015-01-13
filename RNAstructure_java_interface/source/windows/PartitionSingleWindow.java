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
 * A class responsible for initializing and running the partition module, which
 * calculates the partition function for a single strand of nucleic acids.
 *
 * @author Jessica S. Reuter
 */
public class PartitionSingleWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public PartitionSingleWindow( String acid ) {
		super( acid, acid + " Single Strand Partition Function" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Save File" button, try to select an
		// output file, and if one was selected set its name.
		if( command.equals( "Save File" ) ) {
			int index = 2;
			String file = FileChooser.doSave( files.getFile( index ) )
				.getPartition();
			if( !file.equals( "" ) ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Sequence File" button, get a sequence
		// file, initialize a data structure, and create a default output file
		// name.
		else if( command.equals( "Sequence File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = FileChooser.doOpen().getSequence();
			if( file.equals( "" ) ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildPartitionDataStructure( file, isRNA );
			if( !moduleInitialized( result ) ) { return; }

			// Set the sequence file name and the default output file name in
			// the input panel.
			// Then, enable the menus.
			String defaultOut = replaceExtension( file, "pfs" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			menuBar.enableMenus();
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField seq = FileField.createDisabled( "Sequence File" );
		FileField pfs = FileField.createEnabled( "Save File" );
		FilePanel files = new FilePanel( this, seq, pfs );
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

		// Get the data from the window.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String pfs = files.getFile( 2 );
		if( files.isError() ) { return; }

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runPartition( pfs );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// Show the dot plot.
		DrawingWindow drawing = new DrawingWindow( pfs );
		if( drawing.isError() == false ) { drawing.viewWindow(); }
	}

	@Override
	protected RolloverMenu[] setMenus() {

		// Create the temperature menu.
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();

		// Create the menu that handles forced constraints.
		ConstraintsMenu forced = new ConstraintsMenu( backend );
		forced.addGeneralSection();
		forced.addMaxPairingDistanceSection();
		forced.addSHAPESection();
		forced.addShowResetSection();
		forced.addSaveRestoreSection();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, forced };
	}
}
