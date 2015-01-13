/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import java.awt.Dimension;

import RNAstructure_java_interface.source.menus.ConstraintsMenu;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_utilities.FieldPanel.FilePanel;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.FileField;
import RNAstructure_java_utilities.RadioButtonPanel;

/**
 * A class responsible for initializing and running the OligoScreen module,
 * which calculates thermodynamic properties for a list of oligonucleotides.
 *
 * @author Jessica S. Reuter
 */
public class OligoScreenWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The oligomer chemistry.
	 */
	private static String oligomerChemistry;

	/**
	 * Constructor.
	 * <br><br>
	 * This module has DNA hard-coded as its nucleic acid type to start with,
	 * but the nucleic acid type is dynamic and can be changed in the window.
	 */
	public OligoScreenWindow() {
		super( "", "OligoScreen" );
		isRNA = false;
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Oligomer List" button, get a list
		// file, initialize a data structure, and create a default output
		// file name.
		if( command.equals( "Oligomer List" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			String file = FileChooser.doOpen().getList();
			if( file.equals( "" ) ) { return; }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildOligoScreenDataStructure();
			if( !moduleInitialized( result ) ) { return; }

			// Set the list file name and the default output file name in the
			// input panel.
			// Then, enable the menus.
			String defaultOut = replaceExtension( file, "rep" );
			files.setFile( 1, file );
			files.setFile( 2, defaultOut );
			menuBar.enableMenus();
		}

		// If the action comes from the "Report File" button, try to
		// select an output file, and if one was selected set its name.
		else if( command.equals( "Report File" ) ) {
			int index = 2;
			String rep = FileChooser.doSave( files.getFile( index ) )
				.getReport();
			if( !rep.equals( "" ) ) { files.setFile( index, rep ); }
		}

		// If the action comes from the oligomer chemistry group, set the
		// current oligomer chemistry.
		else if( command.equals( "DNA" ) || command.equals( "RNA" ) ) {
			oligomerChemistry = command;
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField list = FileField.createDisabled( "Oligomer List" );
		FileField report = FileField.createEnabled( "Report File" );
		FilePanel files = new FilePanel( this, list, report );
		files.setPanelWidth( 300 );
		files.makePanel();

		// Create the oligo chemistry panel.
		RadioButtonPanel chem = RadioButtonPanel.makeHorizontal(
			"Oligomer Chemistry", this, "DNA", "RNA" );
		int height = chem.getPreferredSize().height;
		chem.setPreferredSize( new Dimension( 200, height ) );
		if( oligomerChemistry != null ) {
			chem.setSelectedButton( oligomerChemistry );
		}

		// Add the components in their proper places.
		setGrid( 2, 1 );
		setFillHorizontal(); 
		placeComponent( 0, 0, files );
		setGrid( 1, 1 );
		setAnchorCenter();
		setInsets( 0, 10, 0, 0 );
		placeComponent( 0, 1, chem );
		makeStartButton( 1, 1 );
	}

	@Override
	protected void runMainCalculation() {

		// Get all input controls from the window.
		FilePanel files = (FilePanel)getInputControl( 1 );
		RadioButtonPanel chem = (RadioButtonPanel)getInputControl( 2 );

		// Get the data from the file input panel.
		// If an error occurred while retrieving data, return.
		String list = files.getFile( 1 );
		String report = files.getFile( 2 );
		if( files.isError() ) { return; }

		// Determine the nucleic acid chemistry type.
		isRNA = chem.getSelectedName().equals( "RNA" );

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarIndeterminate();
		String result =
			backend.runOligoScreen( list, report, isRNA );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// Show a message saying that the calculation has finished.
		dialogHandler.makeMessageDialog( "OligoScreen Complete." );
	}

	@Override
	protected RolloverMenu[] setMenus() {
		ConstraintsMenu temperature = new ConstraintsMenu( backend );
		temperature.buildTemperatureMenu();
		return new ConstraintsMenu[]{ temperature };
	}
}
