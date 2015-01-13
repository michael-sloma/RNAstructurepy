/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.runners;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

import RNAstructure_java_utilities.FileChooser;

/**
 * A class that templates the way in which the RNAstructure Java drawing
 * framework is run independently.
 *
 * @author Jessica S. Reuter
 */
public abstract class RunnerTemplate {

	/**
	 * Constructor.
	 *
	 * @param libraryName   The name of the native library connected to this
	 *                      template.
	 */
	public RunnerTemplate( String libraryName ) {
		try {
			UIManager.setLookAndFeel(
				UIManager.getSystemLookAndFeelClassName() );
			System.loadLibrary( libraryName );
		} catch( Exception e ) { e.getStackTrace(); }
	}

	/**
	 * Add a listener to the dialog spawned by this template that exits the
	 * application when the dialog is closed.
	 *
	 * @param dialog   The dialog to close.
	 */
	public void addClosingListener( JDialog dialog ) {
		dialog.addWindowListener( new WindowAdapter() {
			public void windowClosed( WindowEvent e ) { System.exit( 0 ); }
		});
	}

	/**
	 * Select an image file.
	 *
	 * @return   The image file.
	 */
	protected abstract String selectFile();
}
