/*
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.startup;

import java.awt.SplashScreen;
import java.io.Serializable;

import javax.swing.UIManager;

import RNAstructure_java_utilities.MacChecker;

/**
 * A class that starts the RNAstructure interface.
 *
 * @author Jessica S. Reuter
 */
public class RNAstructure
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * The main method.
	 *
	 * @param args   The command line arguments (ignored).
	 */
	public static void main( String[] args ) {
		try {

			// If the OS is a type of Mac, certain things need to be set in
			// order to make the native Java more "Mac-like."
			MacChecker.checkOS();

			// Set the look and feel for the current OS.
			String lookAndFeel =
				UIManager.getSystemLookAndFeelClassName();
			UIManager.setLookAndFeel( lookAndFeel );

			// Call up the splash screen while the shared library loads.
			SplashScreen splash = SplashScreen.getSplashScreen();
			System.loadLibrary( "RNAstructure_GUI" );
			Thread.sleep( 2000 );
			splash.close();

			// Start the GUI.
			new ApplicationRootFrame();
		} catch( Exception e ) {
			System.err.println( "Error loading RNAstructure." );
			e.printStackTrace();
			System.exit( -1 );
		}
	}
}
