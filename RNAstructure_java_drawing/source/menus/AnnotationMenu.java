/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.menus;

import java.util.ArrayList;

import RNAstructure_java_drawing.source.dialogs.StructureDialog;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.RNAstructureMenu;

/**
 * A class that handles a menu which can annotate a drawn structure.
 *
 * @author Jessica S. Reuter
 */
public class AnnotationMenu
	extends RNAstructureMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * The structure panel connected to this menu.
	 */
	private StructureDialog dialog;

	/**
	 * Constructor.
	 *
	 * @param dialog   The structure dialog connected to this menu.
	 */
	public AnnotationMenu( StructureDialog dialog ) {

		// Create the menu.
		super( "Annotations" );
		addMenuItem(
			"Add Probability Annotation",
			"Add probability color annotation to the structure." );
		addMenuItem(
			"Add SHAPE Annotation",
			"Add SHAPE color annotation to the structure." );
		addMenuItem(
			"Remove Annotation",
			"Remove color annotation from this structure." );

		// Connect the structure panel.
		this.dialog = dialog;
	}

	@Override
	protected void doMenuActions( String command ) {

		// Select a probability data file for annotation, if necessary.
		if( command.contains( "Probability" ) ) {
			String file = FileChooser.doOpen().getPartition();
			if( !file.equals( "" ) ) { dialog.setAnnotation( file ); }
		}

		// Remove annotation, if necessary.
		else if( command.contains( "Remove" ) ) { dialog.clearAnnotation(); }

		// Select a SHAPE data file for annotation, if necessary.
		else if( command.contains( "SHAPE" ) ) {
			String file = FileChooser.doOpen().getSHAPE();
			if( !file.equals( "" ) ) { dialog.setAnnotation( file ); }
		}
	}
}
