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
import RNAstructure_java_utilities.HTMLCheckBox;
import RNAstructure_java_utilities.FieldPanel.FilePanel;
import RNAstructure_java_utilities.FieldPanel.NumberPanel;
import RNAstructure_java_utilities.NumberField.FloatField;
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for initializing and running the bifold module, which
 * folds two strands of nucleic acids into their likely hybrid structure.
 *
 * @author Jessica S. Reuter
 */
public class FoldDoubleWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public FoldDoubleWindow( String acid ) {
		super( acid, acid + " Bimolecular Fold" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, try to select a CT
		// file, and if one was selected set its name.
		if( command.equals( "CT File" ) ) {
			int index = 3;
			String file = FileChooser.doSave( files.getFile( index ) ).getCT();
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
		// Note that all this can only be done if the first sequence file has
		// already been selected.
		else if( command.equals( "Sequence File 2" ) ) {

			// If the first sequence file hasn't been selected yet, show an
			// error and return.
			if( files.getFile( 1 ).equals( "" ) ) {
				String message = "Sequence File 1 hasn't been selected yet.";
				dialogHandler.makeErrorDialog( message );
				return;
			}

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
				backend.buildBifoldDataStructure( seq1, seq2, isRNA );
			if( !moduleInitialized( result ) ) { return; }

			// Build the default output file name and set it in the panel.
			String defaultOut = combineFileNames( seq1, seq2, "ct" );
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
		FileField ct = FileField.createEnabled( "CT File" );
		FilePanel files = new FilePanel( this, seq1, seq2, ct );
		files.setPanelWidth( 400 );
		files.makePanel();

		// Create the box that allows the user to create a save file.
		HTMLCheckBox box =
			HTMLCheckBox.createEmptyBox( "Generate Save File" );

		// Create the parameter panel.
		FloatField energy =
			new FloatField( "Max % Energy Difference", 50, 0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 20, 1 );
		IntegerField window = new IntegerField( "Window Size", 0, 0 );
		NumberPanel params = new NumberPanel( energy, structures, window );
		params.setPanelWidth( 300 );
		params.makePanel();

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setFillCenter(); 
		placeComponent( 0, 1, box );
		setGrid( 1, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 2, params );
		makeStartButton( 1, 2 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		HTMLCheckBox box = (HTMLCheckBox)getInputControl( 2 );
		NumberPanel params = (NumberPanel)getInputControl( 3 );

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		files.getFile( 2 );
		String ct = files.getFile( 3 );
		if( files.isError() ) { return; }

		// Get the data from the parameters panel.
		Float percent = ((FloatField)params.getField( 1 )).getValue();
		Integer structs = ((IntegerField)params.getField( 2 )).getValue();
		Integer window = ((IntegerField)params.getField( 3 )).getValue();

		// Get whether a save file should be written, and whether
		// intramolecular pairs should be forbidden.
		boolean save = box.isSelected();
		boolean noPairs = menuBar.getMenu( 4 ).getItem( 0 ).isSelected();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runBifold( ct, percent, structs, window, save, noPairs );
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

		// Create the menu that determines if unimolecular pairs are allowed.
		ConstraintsMenu pairs = new ConstraintsMenu( backend );
		pairs.buildUnimolecularMenu();

		// Create the maximum loop menu.
		ConstraintsMenu loop = new ConstraintsMenu( backend );
		loop.buildMaxLoopMenu();

		// Return the array of variable menus.
		return new ConstraintsMenu[]{ temperature, pairs, loop };
	}
}
