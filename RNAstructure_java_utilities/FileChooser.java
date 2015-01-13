/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.io.File;
import java.io.Serializable;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * A class that creates a file chooser so the user can select a file to open
 * or save.
 *
 * @author Jessica S. Reuter
 */
public class FileChooser
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * An enum that holds the possible values for the file chooser type.
	 */
	private static enum Types { OPEN, SAVE };

	/**
	 * The default file name for this file chooser, if necessary.
	 */
	private String defaultFile;

	/**
	 * The type of file type selector.
	 */
	private Types selectorType;

	/**
	 * Constructor.
	 *
	 * @param type   The type of file selector.
	 * @param file   The default file in this chooser.
	 */
	private FileChooser( Types type, String file )
		throws IllegalArgumentException {
		selectorType = type;
		defaultFile = file;
	}

	/**
	 * Create a file chooser that selects a file to open.
	 *
	 * @return   The file chooser.
	 */
	public static FileChooser doOpen() {
		return new FileChooser( Types.OPEN, "" );
	}

	/**
	 * Create a file chooser that selects a file to save.
	 *
	 * @return   The file chooser.
	 */
	public static FileChooser doSave() {
		return new FileChooser( Types.SAVE, "" );
	}

	/**
	 * Create a file chooser that selects a file to save, and will have a
	 * default file filled in.
	 *
	 * @param file   The default file name.
	 * @return       The file chooser.
	 */
	public static FileChooser doSave( String file ) {
		return new FileChooser( Types.SAVE, file );
	}

	/**
	 * Select an alignment file.
	 *
	 * @return   The selected alignment file.
	 */
	public String getAlignment() { return select( "Alignment", "ali" ); }

	/**
	 * Select a dot bracket file.
	 *
	 * @return   The selected dot bracket file.
	 */
	public String getBracket() { return select( "Bracket", "bracket" ); }

	/**
	 * Select a constraints file.
	 *
	 * @return   The selected constraints file.
	 */
	public String getConstraints() { return select( "Constraint", "con" ); }

	/**
	 * Select a CT file.
	 *
	 * @return   The selected CT file.
	 */
	public String getCT() { return select( "CT", "ct" ); }

	/**
	 * Select a dot plot file.
	 *
	 * @return   The selected dot plot file.
	 */
	public String getDotPlot() { return select( "Dot Plot", "dp" ); }

	/**
	 * Select a Dynalign save file.
	 *
	 * @return   The selected Dynalign save file.
	 */
	public String getDynalign() { return select( "Dynalign Save", "dsv" ); }

	/**
	 * Select a folding save file.
	 *
	 * @return   The selected folding save file.
	 */
	public String getFolding() { return select( "Folding Save", "sav" ); }

	/**
	 * Select a structure helix file.
	 *
	 * @return   The selected structure helix file.
	 */
	public String getHelix() { return select( "Helix", "txt" ); }

	/**
	 * Select an oligo list file.
	 *
	 * @return   The selected oligo list file.
	 */
	public String getList() { return select( "List", "lis" ); }

	/**
	 * Select a general OUT file.
	 *
	 * @return   The general OUT file.
	 */
	public String getOUT() { return select( "OUT", "out" ); }

	/**
	 * Select a partition function save file.
	 *
	 * @return   The partition function save file.
	 */
	public String getPartition() {
		return select( "Partition Function Save", "pfs" );
	}

	/**
	 * Select a Postscript image file.
	 *
	 * @return   The Postscript image file.
	 */
	public String getPostscript() { return select( "Postscript", "ps" ); }

	/**
	 * Select an oligo report file.
	 *
	 * @return   The oligo report file.
	 */
	public String getReport() { return select( "Report", "rep" ); }

	/**
	 * Select a sequence file.
	 *
	 * @return   The sequence file.
	 */
	public String getSequence() {
		String[] fileTypes = new String[]{ "FASTA", "Sequence" };
		String[] extensions = new String[]{ "fasta", "seq" };
		return select( fileTypes, extensions, false );
	}

	/**
	 * Select a sequence file in its extended context, with more than its
	 * usual possible file types.
	 *
	 * @return   The sequence file.
	 */
	public String getSequenceExtended() {
		String[] fileTypes =
			new String[]{ "FASTA",  "Genbank", "Text", "Sequence"};
		String[] extensions = new String[]{ "fasta","gen", "txt","seq" };
		return select( fileTypes, extensions, true );//added true to allow all file types
	}

	/**
	 * Select a SHAPE file.
	 *
	 * @return   The SHAPE file.
	 */
	public String getSHAPE() { return select( "SHAPE", "shape" ); }

	/**
	 * Select an SVG file.
	 *
	 * @return   The SVG file.
	 */
	public String getSVG() { return select( "SVG", "svg" ); }

	/**
	 * Select a file that has one possible input or output type.
	 *
	 * @param fileType    The type of file to select.
	 * @param extension   The extension of the file types
	 * @return            The name of the selected file.
	 */
	private String select( String fileType, String extension ) {
		return select( new String[]{ fileType }, new String[]{ extension }, false );
	}

	/**
	 * Select a file that has multiple possible input or output types.
	 *
	 * @param fileTypes    The types of file to select.
	 * @param extensions   The extensions of the file types.
	 * @param ShowAll      Indicates whether the "All Files" should be included.  By default, false.
	 * @return             The name of the selected file.
	 */
	private String select( String[] fileTypes, String[] extensions , Boolean ShowAll) {

		// If the number of file types isn't equal to the number of
		// extensions, return the empty string.
		if( fileTypes.length != extensions.length ) { return ""; }

		// Get the current file chooser home directory as a file.
		File currentDirectory =
			new File( System.getProperty( "user.home" ) );

		// Create the file chooser and set its defaults.
		JFileChooser chooser = new JFileChooser();
		chooser.setCurrentDirectory( currentDirectory );
		chooser.setFileSelectionMode( JFileChooser.FILES_ONLY );
		chooser.setAcceptAllFileFilterUsed( ShowAll );
		if( !defaultFile.equals( "" ) ) {
			chooser.setSelectedFile( new File( defaultFile ) );
		}

		// Add the file filter(s).
		int length = fileTypes.length;
		for( int i = 1; i <= length; i++ ) {
			String fullDescription =
				fileTypes[i-1] + " Files (*." + extensions[i-1] + " )";
			chooser.addChoosableFileFilter( new FileNameExtensionFilter(
				fullDescription, extensions[i-1] ) );
		}

		// Show the file dialog, and then if no file was selected, return the
		// empty string.
		int approvalValue = ( selectorType == Types.OPEN ) ?
			chooser.showOpenDialog( null ) : chooser.showSaveDialog( null );
		if( approvalValue != JFileChooser.APPROVE_OPTION ) { return ""; }

		// Get the file, and make sure it ends with the proper extension.
		String file = chooser.getSelectedFile().getAbsolutePath().trim();
		String description = chooser.getFileFilter().getDescription();
		int start = description.lastIndexOf( "." ) + 1;
		int end = description.lastIndexOf( ")" );
		if( ( start != -1 ) && ( end != -1 ) ) {
			String extension = description.substring( start, end ).trim();
			if( !file.endsWith( extension ) ) {
				file = file.concat( "." + extension );
			}
		}

		// Set the file chooser directory based on the file name.
		File fileObject = new File( file );
		String currentDirectoryString = fileObject.getParentFile().getPath();
		System.setProperty( "user.home", currentDirectoryString );

		// Return the name of the file.
		return file;
	}
}
