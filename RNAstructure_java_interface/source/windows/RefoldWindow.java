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
import RNAstructure_java_utilities.FieldPanel.NumberPanel;
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for initializing and running the refold module, which
 * refolds a previous folding calculation from a save file.
 *
 * @author Jessica S. Reuter
 */
public class RefoldWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since the folding save file knows its nucleic acid type, the
	 * exact nucleic acid type this window uses doesn't matter, and the
	 * nucleic acid type is specified as an empty string.
	 */
	public RefoldWindow() {
		super( "", "Refold From Save File" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );
		NumberPanel params = (NumberPanel)getInputControl( 2 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( command.equals( "CT File" ) ) {
			int index = 2;
			String file = FileChooser.doSave( files.getFile( index ) ).getCT();
			if( !file.equals( "" ) ) { files.setFile( index, file ); }
		}

		// If the action comes from the "Save File" button, get a SAV file,
		// initialize a data structure, and create a default output file name.
		else if( command.equals( "Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = FileChooser.doOpen().getFolding();
			if( file.equals( "" ) ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildRefoldSingleDataStructure( file );
			if( !moduleInitialized( result ) ) { return; }

			// Set the folding save file name and the default output file name
			// in the input panel.
			// Then, enable the menus.
			String defaultOut = replaceExtension( file, "ct" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			menuBar.enableMenus();

			// Reset the window size text field based on the sequence length.
			int size = backend.getRefoldWindowSize();
			((IntegerField)params.getField( 3 ))
				.resetField( size, 0, Integer.MAX_VALUE );
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField save = FileField.createDisabled( "Save File" );
		FileField ct = FileField.createEnabled( "CT File" );
		FilePanel files = new FilePanel( this, save, ct );
		files.setPanelWidth( 300 );
		files.makePanel();

		// Create the parameter panel.
		IntegerField energy =
			new IntegerField( "Max % Energy Difference", 10, 1 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 20, 1 );
		IntegerField window = new IntegerField( "Window Size", 0, 0 );
		NumberPanel params = new NumberPanel( energy, structures, window );
		params.setPanelWidth( 250 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setGrid( 1, 1 );
		placeComponent( 0, 1, params );
		makeStartButton( 1, 1 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		NumberPanel params = (NumberPanel)getInputControl( 2 );

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String ct = files.getFile( 2 );
		if( files.isError() ) { return; }

		// Get the data from the parameters panel.
		Integer percent = ((IntegerField)params.getField( 1 )).getValue();
		Integer structures = ((IntegerField)params.getField( 2 )).getValue();
		Integer window = ((IntegerField)params.getField( 3 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runRefold( ct, percent, structures, window );
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
	protected RolloverMenu[] setMenus() {

		// Create the temperature menu.
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();

		// Create the menu that handles forced constraints.
		ConstraintsMenu forced = new ConstraintsMenu( backend );
		forced.addShowResetSection();

		// Create the maximum loop menu.
		ConstraintsMenu loop = new ConstraintsMenu( backend );
		loop.buildMaxLoopMenu();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, forced, loop };
	}
}
