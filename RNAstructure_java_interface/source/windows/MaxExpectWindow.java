/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.FileField;
import RNAstructure_java_utilities.FieldPanel.FilePanel;
import RNAstructure_java_utilities.FieldPanel.NumberPanel;
import RNAstructure_java_utilities.NumberField.DoubleField;
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for initializing and running the MaxExpect module,
 * which generates a group of structures containing pairs with the maximum
 * probability of being accurate.
 *
 * @author Jessica S. Reuter
 */
public class MaxExpectWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public MaxExpectWindow( String acid ) {
		super( acid, acid + " Maximium Expected Accuracy" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( command.equals( "CT File" ) ) {
			int index = 2;
			String file = FileChooser.doSave( files.getFile( index ) ).getCT();
			if( !file.equals( "" ) ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Partition Function Save File" button,
		// get a PFS file, initialize a data structure, and create a default
		// output file name.
		else if( command.equals( "Partition Function Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = FileChooser.doOpen().getPartition();
			if( file.equals( "" ) ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildMaxExpectDataStructure( file );
			if( !moduleInitialized( result ) ) { return; }

			// Set the partition function save file name and the default
			// output file name in the input panel.
			// Then, enable the menus.
			String defaultOut = replaceExtension( file, "ct" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			menuBar.enableMenus();
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField pfs =
			FileField.createDisabled( "Partition Function Save File" );
		FileField ct = FileField.createEnabled( "CT File" );
		FilePanel files = new FilePanel( this, pfs, ct );
		files.setPanelWidth( 400 );
		files.makePanel();

		// Create the gamma weight panel.
		DoubleField gamma = new DoubleField( "Gamma", 1 );
		NumberPanel weight = new NumberPanel( gamma );
		weight.setPanelWidth( 200 );
		weight.makePanel();

		// Create the suboptimal structures parameter panel.
		DoubleField score =
			new DoubleField( "Max % Score Difference", 50, 0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 1000, 1 );
		IntegerField window = new IntegerField( "Window Size", 5, 0 );
		NumberPanel params = new NumberPanel( score, structures, window );
		params.setPanelWidth( 300 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setFillCenter();
		placeComponent( 0, 1, weight );
		setGrid( 1, 1 );
		setFillHorizontal();
		placeComponent( 0, 2, params );
		makeStartButton( 1, 2 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		NumberPanel weight = (NumberPanel)getInputControl( 2 );
		NumberPanel params = (NumberPanel)getInputControl( 3 );

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String ct = files.getFile( 2 );
		if( files.isError() ) { return; }

		// Get the data from the gamma panel and the parameters panel.
		Double gamma = ((DoubleField)weight.getField( 1 )).getValue();
		Double score = ((DoubleField)params.getField( 1 )).getValue();
		Integer structures = ((IntegerField)params.getField( 2 )).getValue();
		Integer window = ((IntegerField)params.getField( 3 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runMaxExpect( ct, score, structures, window, gamma );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// If the user wants to draw structures, draw them.
		String query = "Do you want to draw structures?";
		String doDraw = new SimpleDialogHandler().makeTwoChoicesDialog( query );
		if( doDraw.equals( "OK" ) ) {
			DrawingWindow drawing = new DrawingWindow( ct );
			if( drawing.isError() == false ) { drawing.viewWindow(); }
		}
	}

	@Override
	protected RolloverMenu[] setMenus() { return null; }
}
