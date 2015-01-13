/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.menus;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JLabel;

import RNAstructure_java_drawing.source.dialogs.StructureDialog;
import RNAstructure_java_drawing.source.proxy.StructureBackend;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.RNAstructureMenu;
import RNAstructure_java_utilities.ValueSelectionDialog;

/**
 * A class that handles a menu which can manipulate a drawn structure.
 *
 * @author Jessica S. Reuter
 */
public class StructureMenu
	extends RNAstructureMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * The structure dialog connected to this menu.
	 */
	private StructureDialog dialog;

	/**
	 * Constructor.
	 *
	 * @param dialog   The structure dialog connected to this menu.
	 */
	public StructureMenu( StructureDialog dialog ) {

		// Create the menu.
		super( "Draw" );
		addMenuItem(
			"Go to Structure...",
			"Switch to a selected structure." );
		addMenuItem(
			"Zoom",
			"Zoom this structure." );
		addSeparator();
		JCheckBoxMenuItem clockwise = addCheckBoxMenuItem(
			"Render Clockwise/Counterclockwise",
			"Checked if a structure is rendered clockwise." );
		clockwise.setSelected( StructureDialog.CLOCKWISE_STRUCTURE );
		JCheckBoxMenuItem circle = addCheckBoxMenuItem(
			"Render Nucleotides Circled",
			"Checked if nucleotides are surrounded by circles." );
		circle.setSelected( StructureDialog.CIRCLED_NUCS );
		addSeparator();
		addMenuItem(
			"Write Dot Bracket File",
			"Write a dot bracket file of all structures in this window." );
		addMenuItem(
			"Write Helix (Text) File",
			"Write a text file of the helices in this structure." );
		addMenuItem(
			"Write Postscript File",
			"Write the current structure to a Postscript image file." );
		addMenuItem(
			"Write SVG File",
			"Write the current structure to an SVG image file." );

		// Connect the structure dialog.
		this.dialog = dialog;
	}

	@Override
	protected void doMenuActions( String command ) {

		// If the command comes from the circled nucleotides item, switch
		// circling of nucleotides on or off.
		if( command.endsWith( "Circled" ) ) { dialog.switchCircled(); }

		// If the command is to flip the image, do so.
		else if( command.contains( "Clockwise" ) ) { dialog.flip(); }

		// If the command is to write a particular type of output file, do so.
		else if( command.contains( "File" ) ) {

			// Write a dot bracket file, if necessary.
			if( command.contains( "Dot Bracket" ) ) {
				String file = FileChooser.doSave().getBracket();
				if( !file.equals( "" ) ) { dialog.writeDotBracketFile( file ); }
			}

			// Write a helix file, if necessary.
			else if( command.contains( "Helix" ) ) {
				String file = FileChooser.doSave().getHelix();
				if( !file.equals( "" ) ) { dialog.writeHelixFile( file ); }
			}

			// Write a Postscript file, if necessary.
			else if( command.contains( "Postscript" ) ) {
				String file = FileChooser.doSave().getPostscript();
				if( !file.equals( "" ) ) { dialog.writePostscriptFile( file ); }
			}

			// Write an SVG file, if necessary.
			else if( command.contains( "SVG" ) ) {
				String file = FileChooser.doSave().getSVG();
				if( !file.equals( "" ) ) { dialog.writeSVGFile( file ); }
			}
		}

		// If the command comes from the "Go To Structure..." item, switch
		// structures in the structure window.
		else if( command.equals( "Go to Structure..." ) ) { new MoveDialog(); }

		// If the command is to zoom the image, show a zooming dialog.
		else if( command.equals( "Zoom" ) ) { dialog.viewZoomDialog(); }
	}

	/**
	 * An inner class which creates a dialog that moves to a structure.
	 *
	 * @author Jessica S. Reuter
	 */
	private class MoveDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public MoveDialog() {

			// Call the superclass.
			super(  "Go to Structure" );

			// Get the current and max structure numbers from the dialog label.
			JLabel label = (JLabel)dialog.getContentPane().getComponent( 0 );
			String[] words = label.getText().split( " " );
			Integer current = Integer.parseInt( words[1] );
			Integer max = Integer.parseInt( words[3] );

			// Create the structure field and build the dialog with it.
			IntegerField structureField =
				new IntegerField( "Structure Number", current, 1, max );
			buildDialog( "OK", structureField );
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					dialog.setStructureNumber(
						Integer.parseInt( fields[0].getText() ) );
					dialog.zoomImage();
					dialog.repaint();
					dispose();
				}
			};
		}
	}
}
