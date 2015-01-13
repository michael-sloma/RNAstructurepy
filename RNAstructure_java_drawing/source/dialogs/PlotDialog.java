/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_drawing.source.dialogs;

import java.awt.Graphics2D;
import java.awt.Color;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;

import RNAstructure_java_drawing.source.menus.DotsMenu;
import RNAstructure_java_drawing.source.menus.PlotMenu;
import RNAstructure_java_drawing.source.proxy.DotPlotBackend;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A dialog that holds dot plot images.
 *
 * @author Jessica S. Reuter
 */
public class PlotDialog
	extends ImageDialog {
	private static final long serialVersionUID = 20120802;

	/**
	 * The back end object that reads data for dot plots.
	 */
	private DotPlotBackend backend;

	/**
	 * The strand being drawn.
	 */
	private int strand;

	/**
	 * The length of the plot.
	 */
	private int plotLength;

	/**
	 * File Constructor.
	 * 
	 * @param file   The dot plot data file to draw.
	 */
	public PlotDialog( String file ) {

		// Call the superclass and set the menu bar.
		super( file );
		createDefaultMenuBar(
			new DotsMenu( this ),
			new PlotMenu( this ) );

		// Set the strand to 1, because if no strand is specified, there is
		// only one strand available to begin with.
		strand = 1;

		// Set the top caption to generic text to start with.
		final String caption = "DOT PLOT";
		setTopCaption( caption );

		// Add a mouse listener that allows the user to get information about
		// specific dots. The listener also resets the caption to the default
		// if the area clicked wasn't a dot.
		mainPanel.addMouseListener( new MouseAdapter() {

			@Override
			public void mousePressed( MouseEvent e ) {
				String dotDataString =
					backend.getDotDataByLocation( e.getX(), e.getY(), scale );
				if( dotDataString.equals( "" ) ) { setTopCaption( caption ); }
				else { setTopCaption( dotDataString ); }
			}
		});
	}

	/**
	 * Strand Constructor.
	 * 
	 * @param file     The dot plot data file to draw.
	 * @param strand   The strand from the data file to draw.
	 */
	public PlotDialog( String file, int strand ) {
		this( file );
		this.strand = strand;
	}

	@Override
	protected void createDrawnImage( Graphics2D g2 ) {

		// Draw the grid lines.
		g2.setColor( Color.BLACK );
		int gridLineCounter = 1;
		while( true ) {

			// Get the next grid line.
			// If it doesn't exist, break out of the loop.
			String dataString = backend.getGridLine( gridLineCounter );
			if( dataString.equals( "" ) ) { break; }

			// Split the grid line into data.
			String[] data = dataString.split( " " );
			Double x1 = Double.parseDouble( data[0] );
			Double y1 = Double.parseDouble( data[1] );
			Double x2 = Double.parseDouble( data[2] );
			Double y2 = Double.parseDouble( data[3] );

			// Draw the next grid line.
			g2.draw( new Line2D.Double( x1, y1, x2, y2 ) );

			// Add a label, if necessary.
			if( data.length > 4 ) {

				// Get the label data.
				Integer x = Integer.parseInt( data[5] );
				Integer y = Integer.parseInt( data[6] );
				String text = data[4];

				// Draw the number label.
				g2.drawString( text, x, y );
			}

			// Increment the grid line counter.
			gridLineCounter++;
		}

		// Draw each possible dot.
		for( int i = 1; i <= plotLength; i++ ) {
			for( int j = i; j <= plotLength; j++ ) {

				// Get the next dot's data.
				String data = backend.getDotData( i, j );
				if( !data.equals( "" ) ) {

					// Split the data into its pieces.
					String[] dataArray = data.split( " " );
					Integer x = Integer.parseInt( dataArray[0] );
					Integer y = Integer.parseInt( dataArray[1] );
					Integer red = Integer.parseInt( dataArray[2] );
					Integer green = Integer.parseInt( dataArray[3] );
					Integer blue = Integer.parseInt( dataArray[4] );

					// Draw the dot in its appropriate color.
					g2.setColor( new Color( red, green, blue ) );
					g2.fillRect( x, y, 3, 3 );
				}
			}
		}
	}

	/**
	 * Get the range of colors allowed in this dialog.
	 *
	 * @return   The color range.
	 */
	public Integer[] getColors() {
		String[] colorRange = backend.getColors().split( " " );
		Integer current = Integer.parseInt( colorRange[0] );
		Integer min = Integer.parseInt( colorRange[1] );
		Integer max = Integer.parseInt( colorRange[2] );
		return new Integer[]{ current, min, max };
	}

	/**
	 * Get the current bounds allowed in this dialog.
	 *
	 * @return   The current bounds.
	 */
	public Double[] getCurrentBounds() {
		Double min = backend.getCurrentMin();
		Double max = backend.getCurrentMax();
		return new Double[]{ min, max };
	}

	/**
	 * Get the default bounds allowed in this dialog.
	 *
	 * @return   The default bounds.
	 */
	public Double[] getDefaultBounds() {
		Double min = backend.getDefaultMin();
		Double max = backend.getDefaultMax();
		return new Double[]{ min, max };
	}

	@Override
	protected void readImageData() {

		// Read the dot plot data.
		// If an error occurred, set the error and return.
		backend = new DotPlotBackend();
		boolean correctRead = backend.readData( file, strand );
		if( correctRead == false ) {
			error = "Error reading dot plot data.";
			return;
		}

		// Set the plot bounds.
		String[] bounds = backend.getBounds().split( " " );
		maxBoundX = Double.parseDouble( bounds[0] );
		maxBoundY = Double.parseDouble( bounds[1] );

		// Determine the plot length.
		plotLength = backend.getLength();

		// Build the legend.
		buildLegend( backend.getLegendData().split( "\n" ) );

		// Zoom the image.
		zoomImage();
	}

	/**
	 * Set the number of colors in the plot on this dialog.
	 *
	 * @param colors   The number of colors to set.
	 */
	public void setColors( int colors ) {
		backend.setColors( colors );
		buildLegend( backend.getLegendData().split( "\n" ) );
	}

	/**
	 * Set the range of the plot on this dialog.
	 *
	 * @param minimum   The lower bound of the range.
	 * @param maximum   The upper bound of the range.
	 */
	public void setRange( double minimum, double maximum ) {
		backend.setRange( minimum, maximum );
		buildLegend( backend.getLegendData().split( "\n" ) );
	}

	/**
	 * Write a Postscript file from the plot on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writePostscriptFile( String outFile ) {
		backend.writePostscriptFile( outFile );
		String message = "Postscript file written.";
		new SimpleDialogHandler().makeMessageDialog( message );
	}

	/**
	 * Write a probable structures file from the plot on this dialog.
	 *
	 * @param outFile          The file to write.
	 * @param showStructures   True if structures should be shown, false if not.
	 * @return                 True if structures were written correctly, false
	 *                         if not.
	 */
	public boolean writeStructuresFile( String outFile, boolean showStructures ) {

		// Write the structures file and create a new dialog handler.
		String result = backend.writeStructuresFile( file, outFile );
		SimpleDialogHandler handler = new SimpleDialogHandler();

		// Show an error and return if a problem occurred writing the
		// structures file.
		if( !result.equals( "" ) ) {
			handler.makeErrorDialog(
				"Error writing probable structures file." );
			return false;
		}

		// If the option should be given to show structures, give the user a
		// chance to draw them immediately, and show structures if the user
		// requests them.
		if( showStructures == true ) {
			String query = "Do you want to draw structures?";
			String doDraw = handler.makeTwoChoicesDialog( query );
			if( doDraw.equals( "OK" ) ) {
				StructureDialog dialog = new StructureDialog( outFile );
				if( dialog.isError() == false ) { dialog.viewDialog(); }
			}
		}

		// If the method got here, everything went well, so return true.
		return true;
	}

	/**
	 * Write an SVG file from the plot on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeSVGFile( String outFile ) {
		backend.writeSVGFile( outFile );
		new SimpleDialogHandler().makeMessageDialog( "SVG file written." );
	}

	/**
	 * Write a text file from the plot on this dialog.
	 *
	 * @param outFile   The file to write.
	 */
	public void writeTextFile( String outFile ) {
		backend.writeTextFile( outFile );
		String message = "Dot plot text file written.";
		new SimpleDialogHandler().makeMessageDialog( message );
	}
}
