/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.menus;

import java.awt.Desktop;
import java.net.URI;

import RNAstructure_java_utilities.ImageGrabber;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class that creates a "Help" menu.
 *
 * @author Jessica S. Reuter
 */
public class HelpMenu
	extends RolloverMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 */
	public HelpMenu() {
		super( "Help" );
		addMenuItem( "Help Topics", "Get online help.", "F1" );
		addSeparator();
		addMenuItem( "About RNAstructure...",
			"Display program information, version number, and copyright." );
	}

	@Override
	protected void doMenuActions( String command ) {
		SimpleDialogHandler dialog = new SimpleDialogHandler();

		// If the command is to show the about screen, do so.
		if( command.startsWith( "About" ) ) {
			try {
				dialog.makeMessageDialog(
					ImageGrabber.getImageLabel( "/images/Splash.gif" ) );
			} catch( Exception e ) {
				dialog.makeErrorDialog( "Error showing about screen." );
			}
		}

		// Otherwise, show the help in a browser window.
		else {
			try {
				String page =
					"http://rna.urmc.rochester.edu/GUI/html/Contents.html";
				Desktop.getDesktop().browse( new URI( page ) );
			} catch( Exception e ) {
				dialog.makeErrorDialog( "Error showing online help." );
			}
		}
	}
}
