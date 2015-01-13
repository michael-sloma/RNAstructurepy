/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.awt.Color;

import javax.swing.JTextField;

/**
 * A class that creates a text field which holds a file name.
 *
 * @author Jessica S. Reuter
 */
public class FileField
	extends JTextField {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 * 
	 * @param text      The name of the field.
	 * @param enabled   Whether the field is enabled.
	 */
	private FileField( String text, boolean enabled ) {
		setName( text );
		if( !enabled ) {
			setEditable( false );
			setBackground( Color.WHITE );
		}
	}

	/**
	 * Create a disabled file input field.
	 *
	 * @return   The field.
	 */
	public static FileField createDisabled( String text ) {
		return new FileField( text, false );
	}

	/**
	 * Create an enabled file input field.
	 *
	 * @return   The field.
	 */
	public static FileField createEnabled( String text ) {
		return new FileField( text, true );
	}

	/**
	 * Get the file name in the field.
	 *
	 * @return   The file name.
	 */
	public String getFile() { return getText().trim(); }
}
