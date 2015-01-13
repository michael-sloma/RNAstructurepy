/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.runners;

import java.io.Serializable;

import javax.swing.JOptionPane;

import RNAstructure_java_drawing.source.dialogs.StructureDialog;
import RNAstructure_java_utilities.FileChooser;

/**
 * A class that runs RNAstructure Java structure drawing independently of the
 * main RNAstructure application.
 *
 * @author Jessica S. Reuter
 */
public class StructureRunner
	extends RunnerTemplate implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 */
	public StructureRunner() { super( "DrawingStructures" ); }

	@Override
	protected String selectFile() {
		JOptionPane.showMessageDialog( null,
			"Please select a structure file to draw using the " +
			"following file chooser." );
		return FileChooser.doOpen().getCT();
	}

	/**
	 * The main method.
	 *
	 * @param args   The command line arguments (ignored).
	 */
	public static void main( String args[] ) {
		StructureRunner template = new StructureRunner();
		String file = template.selectFile();
		if( !file.equals( "" ) ) {
			StructureDialog dialog = new StructureDialog( file );
			template.addClosingListener( dialog );
			dialog.viewDialog();
		}
	}
}
