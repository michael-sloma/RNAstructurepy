/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.menus;

import RNAstructure_java_drawing.source.dialogs.PlotDialog;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.RNAstructureMenu;

/**
 * A class that handles a menu which can manipulate a drawn dot plot.
 *
 * @author Jessica S. Reuter
 */
public class PlotMenu
	extends RNAstructureMenu {
	private static final long serialVersionUID = 20120615;

	/**
	 * The plot dialog connected to this menu.
	 */
	private PlotDialog dialog;

	/**
	 * Constructor.
	 *
	 * @param dialog   The structure dialog connected to this menu.
	 */
	public PlotMenu( PlotDialog dialog ) {

		// Create the menu.
		super( "Output Plot" );
		addMenuItem(
			"Write Dot Plot File",
			"Write visible dots to a text file." );
		addMenuItem(
			"Write Postscript File",
			"Write the visible plot to a Postscript image file." );
		addMenuItem(
			"Write SVG File",
			"Write the visible plot to an SVG image file." );
		if( dialog.getFile().endsWith( "pfs" ) ) {
			addSeparator();
			addMenuItem(
				"Write Probable Structures File",
				"Write a CT file containing structures which are composed of " +
					"probable pairs of different levels." );
		}

		// Connect the plot dialog.
		this.dialog = dialog;
	}

	@Override
	protected void doMenuActions( String command ) {

		// Write a dot plot file, if necessary.
		if( command.contains( "Dot Plot" ) ) {
			String file = FileChooser.doSave().getDotPlot();
			if( !file.equals( "" ) ) { dialog.writeTextFile( file ); }
		}

		// Write a Postscript file, if necessary.
		else if( command.contains( "Postscript" ) ) {
			String file = FileChooser.doSave().getPostscript();
			if( !file.equals( "" ) ) { dialog.writePostscriptFile( file ); }
		}

		// Write a probable structures file, if necessary.
		else if( command.equals( "Write Probable Structures File" ) ) {
			String file = FileChooser.doSave().getCT();
			if( !file.equals( "" ) ) {
				dialog.writeStructuresFile( file, true );
			}
		}

		// Write an SVG file, if necessary.
		else if( command.contains( "SVG" ) ) {
			String file = FileChooser.doSave().getSVG();
			if( !file.equals( "" ) ) { dialog.writeSVGFile( file ); }
		}
	}
}
