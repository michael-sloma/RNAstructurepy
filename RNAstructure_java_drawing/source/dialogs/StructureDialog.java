/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.dialogs;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JLabel;
import javax.swing.JMenuBar;

import RNAstructure_java_drawing.source.menus.AnnotationMenu;
import RNAstructure_java_drawing.source.menus.StructureMenu;
import RNAstructure_java_drawing.source.proxy.StructureBackend;
import RNAstructure_java_utilities.MacChecker;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A dialog that holds structure images.
 *
 * @author Jessica S. Reuter
 */
public class StructureDialog
	extends ImageDialog {
	private static final long serialVersionUID = 20120802;

	/**
	 * Boolean that tells whether a structure's nucleotides are circled.
	 */
	public static boolean CIRCLED_NUCS = true;

	/**
	 * Boolean that tells whether a structure is rendered clockwise.
	 */
	public static boolean CLOCKWISE_STRUCTURE = true;

	/**
	 * The back end object that reads data for structures.
	 */
	private StructureBackend backend;

	/**
	 * The array of data for a particular structure.
	 */
	private String[] data;

	/**
	 * Whether the current structure has pairs.
	 */
	private boolean hasPairs = false;

	/**
	 * The structure number viewed on this panel.
	 */
	private int number;

	/**
	 * The array of regexes used to pull data about a structure.
	 */
	private String[] regex;

	/**
	 * Constructor.
	 * 
	 * @param file   The structure file to draw.
	 */
	public StructureDialog( String file ) {

		// Call the superclass and set the menu bar.
		super( file );
		createDefaultMenuBar(
			new AnnotationMenu( this ),
			new StructureMenu( this ) );

		// Create the pattern strings.
		regex = new String[6];
		regex[0] =
			"Placed at\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[1] =
			"Paired to\\: " +
			"([0-9]+)";
		regex[2] =
			"Control point\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[3] =
			"Backbone stretches to: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[4] =
			"Label placed at: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\)";
		regex[5] =
			"Character\\: ([A-Za-z]).*" +
			"Placed at\\: " +
			"\\(([0-9]+\\.*[0-9]*)\\,([0-9]+\\.*[0-9]*)\\).*" +
			"Color \\(RGB\\)\\: " +
			"([0-9]+)\\,([0-9]+)\\,([0-9]+)";

		// Add a key listener that allows the user to move between structures.
		addKeyListener( new KeyAdapter() {
			private boolean isControl = false;

			@Override
			public void keyPressed( KeyEvent e ) {
				isControl = e.isControlDown();
				if( MacChecker.isMac() ) {
					if( isControl == false ) { isControl = e.isMetaDown(); }
				}
			}

			@Override
			public void keyReleased( KeyEvent e ) {
				int keyCode = e.getKeyCode();
				boolean enter =
					( keyCode == KeyEvent.VK_UP ) ||
					( keyCode == KeyEvent.VK_DOWN );

				if( enter && isControl ) {
					
					//Current structure is being taken from structure label!
					JLabel label = (JLabel)getContentPane().getComponent( 0 );
					String[] words = label.getText().split( " " );
					Integer current = Integer.parseInt( words[1] );
					Integer max = Integer.parseInt( words[3] );

					if( keyCode == KeyEvent.VK_UP ) { current++; }
					else { current--; }

					if( current < 1 ) { current = 1; }
					else if( current > max ) { current = max; }
					setStructureNumber( current );
				}

				isControl = false;
			}
		});

		// Set the structure number to the first.
		setStructureNumber( 1 );
	}

	/**
	 * Clear annotation from this dialog.
	 */
	public void clearAnnotation() {
		backend.removeAnnotation();
		buildLegend();
		refresh();
	}

	@Override
	protected void createDrawnImage( Graphics2D g2 ) {

		// Set the drawing color to black.
		g2.setColor( Color.black );

		// If the structure doesn't have pairs, paint a message stating this
		// and return.
		if( hasPairs == false ) {
			g2.setFont( g2.getFont().deriveFont( Font.BOLD, 20 ) );
			g2.drawString( "This structure contains no pairs.", 30, 60 );
			return;
		}

		// Set thick lines for the pairs.
		g2.setStroke( new BasicStroke( (float)3 ) );

		// Draw the pairs, where necessary.
		for( int i = 1; i < data.length; i++ ) {

			// Determine whether the next nucleotide is paired.
			// If it is, get the pair data.
			Matcher m = Pattern.compile( regex[1] ).matcher( data[i] );
			int pair = ( m.find() ) ? Integer.parseInt( m.group( 1 ) ) : 0;
			if( pair != 0 ) {

				// Use matchers to find the coordinates and control point of
				// the pair.
				Matcher m1 = Pattern.compile( regex[0] ).matcher( data[i] );
				m1.find();
				Matcher m2 = Pattern.compile( regex[2] ).matcher( data[i] );
				m2.find();
				Matcher m3 = Pattern.compile( regex[0] ).matcher( data[pair] );
				m3.find();

				// Convert the coordinate values to doubles.
				Double x1 = Double.parseDouble( m1.group( 1 ) );
				Double y1 = Double.parseDouble( m1.group( 2 ) );
				Double centerX = Double.parseDouble( m2.group( 1 ) );
				Double centerY = Double.parseDouble( m2.group( 2 ) );
				Double x2 = Double.parseDouble( m3.group( 1 ) );
				Double y2 = Double.parseDouble( m3.group( 2 ) );

				// Draw the pair as a path.
				Path2D arc = new Path2D.Double();
				arc.moveTo( x1, y1 );
				arc.curveTo( x1, y1, centerX, centerY, x2, y2 );
				g2.draw( arc );
			}
		}

		// Set lines of normal thickness after pairs are done being drawn.
		g2.setStroke( new BasicStroke( (float)1 ) );

		// Draw backbone and nucleotide label lines, where necessary.
		for( int i = 1; i < data.length; i++ ) {

			// Check to make sure coordinates exist for the next data point.
			// If they do exist, draw lines.
			Matcher m1 = Pattern.compile( regex[0] ).matcher( data[i] );
			if( m1.find() ) {

				// Get the next nucleotide's location coordinates, as doubles,
				// to use as an anchor for potential backbone or label lines.
				Double nucX = Double.parseDouble( m1.group( 1 ) );
				Double nucY = Double.parseDouble( m1.group( 2 ) );

				// Use matchers to search for coordinates for backbone or
				// label lines.
				Matcher m2 = Pattern.compile( regex[3] ).matcher( data[i] );
				Matcher m3 = Pattern.compile( regex[4] ).matcher( data[i] );

				// Draw a backbone line, if one should exist.
				if( m2.find() ) {
					Double backX = Double.parseDouble( m2.group( 1 ) );
					Double backY = Double.parseDouble( m2.group( 2 ) );
					g2.draw( new Line2D.Double( nucX, nucY, backX, backY ) );
				}

				// Draw a label line, if one should exist.
				if( m3.find() ) {
					Double labelX = Double.parseDouble( m3.group( 1 ) );
					Double labelY = Double.parseDouble( m3.group( 2 ) );
					g2.draw( new Line2D.Double( nucX, nucY, labelX, labelY ) );
				}
			}
		}

		// Set the font bold and big to draw labels.
		g2.setFont( g2.getFont().deriveFont( Font.BOLD, 28 ) );

		// Draw the nucleotide number labels, where they exist.
		for( int i = 1; i < data.length; i++ ) {
			Matcher m = Pattern.compile( regex[4] ).matcher( data[i] );
			if( m.find() ) {

				// Get the label data.
				Integer x =
					new Double( Double.parseDouble( m.group( 1 ) ) ).intValue();
				Integer y =
					new Double( Double.parseDouble( m.group( 2 ) ) ).intValue();
				String text = new Integer( i ).toString();

				// Get the bounds of the string to draw, and the visual bounds
				// of the text itself.
				Rectangle stringBounds =
					g2.getFontMetrics().getStringBounds( text, g2 ).getBounds();
				Rectangle visualBounds =
					g2.getFont()
					.createGlyphVector( g2.getFontRenderContext(), text )
					.getVisualBounds().getBounds();

				// Draw the number label.
				g2.setColor( Color.white );
				g2.fillRect(
					x - ( stringBounds.width / 2 ) - 1,
					y - ( visualBounds.height / 2 ) - 1,
					stringBounds.width + 2,
					visualBounds.height + 2 );
				g2.setColor( Color.black );
				g2.drawString( text,
					x - ( stringBounds.width / 2 ),
					y - ( visualBounds.height / 2 ) - visualBounds.y );
			}
		}

		// Draw the nucleotide letter labels.
		int nucCharRadius = 30;
		int halfCharRadius = nucCharRadius / 2;
		for( int i = 1; i < data.length; i++ ) {

			// Check to make sure label data exists for the next data point.
			// If it does exist, draw the label.
			Matcher m =
				Pattern.compile( regex[5], Pattern.DOTALL ).matcher( data[i] );
			if( m.find() ) {

				// Get the label data.
				Integer x =
					new Double( Double.parseDouble( m.group( 2 ) ) ).intValue();
				Integer y =
					new Double( Double.parseDouble( m.group( 3 ) ) ).intValue();
				String text = m.group( 1 );

				// Get the bounds of the string to draw, and the visual bounds
				// of the text itself.
				Rectangle stringBounds =
					g2.getFontMetrics().getStringBounds( text, g2 ).getBounds();
				Rectangle visualBounds =
					g2.getFont()
					.createGlyphVector( g2.getFontRenderContext(), text )
					.getVisualBounds().getBounds();

				// Determine the color of the nucleotide.
				Integer red = Integer.parseInt( m.group( 4 ) );
				Integer green = Integer.parseInt( m.group( 5 ) );
				Integer blue = Integer.parseInt( m.group( 6 ) );

				// Draw the letter label.
				g2.setColor( Color.white );
				g2.fillOval(
					x - halfCharRadius, y - halfCharRadius,
					nucCharRadius, nucCharRadius );
				if( CIRCLED_NUCS ) {
					g2.setColor( Color.black );
					g2.drawOval(
						x - halfCharRadius, y - halfCharRadius,
						nucCharRadius, nucCharRadius );
				}
				g2.setColor( new Color( red, green, blue ) );
				g2.drawString( text,
					x - ( stringBounds.width / 2 ),
					y - ( visualBounds.height / 2 ) - visualBounds.y );
			}
		}
	}

	/**
	 * Flip the structure image horizontally.
	 */
	public void flip() {
		boolean clockwise = CLOCKWISE_STRUCTURE;
		CLOCKWISE_STRUCTURE = !clockwise;
		backend.flip();
		repaint();
	}

	@Override
	protected void readImageData() {
		backend = new StructureBackend();
		if( ( backend.readStructureData( file ) ) == false ) {
			error = "Error reading structure data.";
			return;
		}
	}

	/**
	 * Set the annotation for this structure drawing.
	 *
	 * @param file   The data file to get structure information from.
	 */
	public void setAnnotation( String file ) {

		// Read the appropriate annotation.
		boolean correctRead =
			( file.endsWith( "pfs" ) ) ?
				backend.addAnnotationProbability( file ) :
			( file.endsWith( "shape" ) ) ?
				backend.addAnnotationSHAPE( file ) :
			false;

		// If annotation was not read correctly, show an error and return.
		if( correctRead == false ) {
			String error = "Error adding annotation.";
			new SimpleDialogHandler().makeErrorDialog( error );
			return;
		}

		// Build the legend.
		String data = backend.getStructureData( 1 );
		data = data.substring( data.indexOf( "Legend" ) ).trim();
		data = data.substring( data.indexOf( ":" ) + 2 );
		data = data.replaceAll( "rgb\\(", "" )
			.replaceAll( "\\)", "" )
			.replaceAll( ",", " " )
			.replaceAll( "\"", "" )
			.replaceAll( " -- ", " " );
		String[] legend = data.split( "\\n" );

		// Place the legend and repaint the dialog.
		buildLegend( legend );
		
		
		//DHM ADDED to FIX BUG where annotations don't appear.
		//Note that keys do not appear until windows is moved or structure number is changed...
		
		//Note how current is coming from the structure label!
		JLabel label = (JLabel)getContentPane().getComponent( 0 );
		String[] words = label.getText().split( " " );
		Integer current = Integer.parseInt( words[1] );
		setStructureNumber(current);
		
		//End added by DHM
		
		refresh();
	}

	/**
	 * Set the structure number.
	 *
	 * @param number   The structure number.
	 */
	public void setStructureNumber( int number ) {

		// Set the structure number and get its data.
		// If an error occurred, show a general error and return.
		this.number = number;
		String dataString = backend.getStructureData( number );
		if( dataString.equals( "" ) ) {
			String message = "Error retrieving structure data.";
			new SimpleDialogHandler().makeErrorDialog( message );
			return;
		}
		data = null;
		data = dataString.split( "Nucleotide" );

		// Set the bounds.
		String bounds =
			dataString.substring( dataString.indexOf( "Max Bounds" ) );
		int xBound1 = bounds.indexOf( "(" ) + 1;
		int xBound2 = bounds.indexOf( "," );
		int yBound1 = xBound2 + 1;
		int yBound2 = bounds.indexOf( ")" );
		maxBoundX = Double.parseDouble( bounds.substring( xBound1, xBound2 ) );
		maxBoundY = Double.parseDouble( bounds.substring( yBound1, yBound2 ) );

		// Set the label text.
		String structureID =
			dataString.substring( 0, dataString.indexOf( "Nuc" ) ).trim();
		String description =
			dataString.substring( dataString.indexOf( "Description" ) );
		description = description.substring( description.indexOf( ":" ) + 2 );
		if( description.contains( "Legend:" ) ) {
			int legendIndex = description.indexOf( "Legend:" );
			description = description.substring( 0, legendIndex );
		}
		setTopCaption( structureID + " -- " + description.trim() );

		// Determine if the structure has pairs.
		// If not, set max bounds equal to the window size and scale to 50%.
		hasPairs = Pattern.compile( regex[1] ).matcher( dataString ).find();
		if( hasPairs == false ) {
			setScale( 50 );
			maxBoundX = preferredDialogWidth;
			maxBoundY = preferredDialogWidth;
		}

		// Calculate scale so the whole image fits in the window at the
		// first glance, and zoom.
		double scaleWidth = preferredDialogWidth - 5;
		double xScale = scaleWidth / maxBoundX;
		double yScale = scaleWidth / maxBoundY;
		scale = Math.min( xScale, yScale );
		zoomImage();
	}

	/**
	 * Switch circled nucleotides on or off.
	 */
	public void switchCircled() {
		boolean circled = CIRCLED_NUCS;
		CIRCLED_NUCS = !circled;
		backend.setNucleotidesCircled( CIRCLED_NUCS );
		refresh();
	}

	/**
	 * Write a dot bracket file of all structures on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeDotBracketFile( String outFile ) {
		SimpleDialogHandler handler = new SimpleDialogHandler();
		String result = backend.writeDotBracketFile( file, outFile );
		if( result.equals( "" ) ) {
			handler.makeMessageDialog( "Dot bracket file written." );
		} else {
			handler.makeErrorDialog( "Error writing dot bracket file." );
		}
	}

	/**
	 * Write a helix file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeHelixFile( String outFile ) {
		SimpleDialogHandler handler = new SimpleDialogHandler();
		String result = backend.writeHelixFile( file, outFile, number );
		if( result.equals( "" ) ) {
			handler.makeMessageDialog( "Helix file written." );
		} else {
			handler.makeErrorDialog( "Error writing helix file." );
		}
	}

	/**
	 * Write a Postscript file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writePostscriptFile( String outFile ) {
		backend.writePostscriptFile( outFile, number );
		new SimpleDialogHandler().makeMessageDialog(
			"Postscript file written." );
	}

	/**
	 * Write an SVG file from the current structure on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeSVGFile( String outFile ) {
		backend.writeSVGFile( outFile, number );
		new SimpleDialogHandler().makeMessageDialog( "SVG file written." );
	}
}
