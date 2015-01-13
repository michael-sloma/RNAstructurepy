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
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for initializing and running the Dynalign refold
 * module, which refolds a previous Dynalign calculation from a save file.
 *
 * @author Jessica S. Reuter
 */
public class DynalignRefoldWindow
	extends ModuleWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 * <br><br>
	 * Note that since the Dynalign save file knows its nucleic acid type, the
	 * exact nucleic acid type this window uses doesn't matter, and the
	 * nucleic acid type is specified as an empty string.
	 */
	public DynalignRefoldWindow() {
		super( "", "Refold From Dynalign Save File" );
	}

	@Override
	protected void executeModuleAction( String command ) {

		// Get any input controls from the window that handle actions.
		FilePanel files = (FilePanel)getInputControl( 1 );

		// If the action comes from the "Alignment File" button, try to select
		// an alignment file, and if one was selected set its name.
		if( command.equals( "Alignment File" ) ) {
			String file = FileChooser.doSave().getAlignment();
			if( !file.equals( "" ) ) { files.setFile( 4, file ); }
		}

		// If the action comes from one of the "CT File" buttons, try to
		// select a CT file, and if one was selected set its name, based on
		// its index.
		else if( command.startsWith( "CT File" ) ) {
			String file = FileChooser.doSave().getCT();
			if( !file.equals( "" ) ) {
				int index = ( command.endsWith( "1" ) ) ? 2 : 3;
				files.setFile( index, file );
			}
		}

		// If the action comes from the "Save File" button, get a Dynalign
		// save file and initialize a data structure.
		else if( command.equals( "Save File" ) ) {

			// Attempt to select the file.
			// If no file was selected, return.
			// Otherwise, set the file in the input panel.
			String file = FileChooser.doOpen().getDynalign();
			if( file.equals( "" ) ) { return; }
			else { files.setFile( 1, file ); }

			// Create a data structure.
			// If an error occurred creating the data structure, show an error
			// and return.
			String result =
				backend.buildRefoldDynalignDataStructure( file );
			if( !moduleInitialized( result ) ) { return; }

			// Enable the menus.
			menuBar.enableMenus();
		}
	}

	@Override
	protected void makeInputControls() {

		// Create the file input panel.
		FileField save = FileField.createDisabled( "Save File" );
		FileField ct1 = FileField.createEnabled( "CT File 1" );
		FileField ct2 = FileField.createEnabled( "CT File 2" );
		FileField align = FileField.createEnabled( "Alignment File" );
		FilePanel files = new FilePanel( this, save, ct1, ct2, align );
		files.setPanelWidth( 400 );
		files.makePanel();

		IntegerField energy =
			new IntegerField( "Max % Energy Difference", 20, 0 );
		IntegerField structures =
			new IntegerField( "Max Number of Structures", 750, 1 );
		IntegerField windowStruct =
			new IntegerField( "Structure Window Size", 0, 0 );
		IntegerField windowAlign =
			new IntegerField( "Alignment Window Size", 0, 0 );
		NumberPanel params =
			new NumberPanel( energy, structures, windowStruct, windowAlign );
		params.setPanelWidth( 300 );
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
		String ct1 = files.getFile( 2 );
		String ct2 = files.getFile( 3 );
		String align = files.getFile( 4 );
		if( files.isError() ) { return; }

		// Get the data from the parameters panel.
		Integer percent = ((IntegerField)params.getField( 1 )).getValue();
		Integer structures = ((IntegerField)params.getField( 2 )).getValue();
		Integer windowStr = ((IntegerField)params.getField( 3 )).getValue();
		Integer windowAli = ((IntegerField)params.getField( 4 )).getValue();

		// Run the calculation.
		// If an error occurred during calculation, return.
		showProgressBarDeterminate();
		String result =
			backend.runDynalignRefold(
				ct1, ct2, align, percent, structures, windowStr, windowAli );
		finishModule( result );
		if( !result.equals( "" ) ) { return; }

		// If the user wants to draw structures, draw them.
		String query = "Do you want to draw structures?";
		String doDraw = new SimpleDialogHandler().makeTwoChoicesDialog( query );
		if( doDraw.equals( "OK" ) ) {
			DrawingWindow one = new DrawingWindow( ct1 );
			DrawingWindow two = new DrawingWindow( ct2 );
			if( ( one.isError() == false ) && ( two.isError() == false ) ) {
				one.viewWindow();
				two.viewWindow();
			}
		}
	}

	@Override
	protected RolloverMenu[] setMenus() { return null; }
}
