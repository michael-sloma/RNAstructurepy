/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import RNAstructure_java_interface.source.menus.ConstraintsMenu;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_utilities.SimpleDialogHandler;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.FileField;
import RNAstructure_java_utilities.HTMLCheckBox;
import RNAstructure_java_utilities.FieldPanel.FilePanel;

/**
 * A class responsible for initializing and running the efn2 (Energy Function
 * 2) module, which calculates the folding free energy of structures.
 *
 * @author Jessica S. Reuter
 */
public class Efn2Window
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param acid   The nucleic acid type.
	 */
	public Efn2Window( String acid ) {
		super( acid, acid + " Efn2" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "CT File" button, get a CT file,
		// initialize a data structure, and create a default output file name.
		if( command.equals( "CT File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = FileChooser.doOpen().getCT();
			if( file.equals( "" ) ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildEfn2DataStructure( file, isRNA );
			if( !moduleInitialized( result ) ) { return; }

			// Set the CT file name and the default output file name in the
			// input panel.
			// Then, enable the menus.
			String defaultOut = replaceExtension( file, "out" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			menuBar.enableMenus();
		}

		// If the action comes from the "Output File" button, try to select an
		// output file, and if one was selected set its name.
		else if( command.equals( "Output File" ) ) {
			int index = 2;
			String file = FileChooser.doSave( files.getFile( index ) ).getOUT();
			if( !file.equals( "" ) ) { files.setFile( index, file ); }
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField ct = FileField.createDisabled( "CT File" );
		FileField out = FileField.createEnabled( "Output File" );
		FilePanel files = new FilePanel( this, ct, out );
		files.setPanelWidth( 300 );
		files.makePanel();

		// Make the thermodynamic details box.
		HTMLCheckBox box =
			HTMLCheckBox.createEmptyBox( "Write Thermodynamic Details File" );

		// Add the components in their proper places.
		setGrid( 1, 1 );
		setFillHorizontal();
		placeComponent( 0, 0, files );
		setGrid( 2, 1 );
		setFillCenter();
		placeComponent( 0, 1, box );
		setAnchorNorth();
		setGrid( 1, 1 );
		makeStartButton( 2, 0 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		HTMLCheckBox box = (HTMLCheckBox)getInputControl( 2 );

		// Get all data from the window.
		// If an error occurred while retrieving data, return.
		files.getFile( 1 );
		String out = files.getFile( 2 );
		if( files.isError() ) { return; }
		boolean writeDetails = box.isSelected();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarIndeterminate();
		String result =
			backend.runEfn2( out, writeDetails );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// Show a message saying that the calculation has finished.
		new SimpleDialogHandler().makeMessageDialog( "Efn2 Complete." );
	}

	@Override
	protected RolloverMenu[] setMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
