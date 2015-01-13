/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import RNAstructure_java_utilities.BorderBuilder;

/**
 * A class responsible for creating and displaying a window that enables a
 * user to input data for one of the RNAstructure modules.
 *
 * @author Jessica S. Reuter
 */
public abstract class ModuleWindow
	extends InternalWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The GridBagConstraints for this window.
	 */
	private GridBagConstraints constraints;

	/**
	 * The layout for this window.
	 */
	private GridBagLayout layout;

	/**
	 * A boolean, true if the calculation handles RNA, false if not.
	 */
	protected boolean isRNA;

	/**
	 * The progress bar dialog for this window.
	 */
	protected final JDialog progressBar;

	/**
	 * Constructor.
	 *
	 * @param acid    The nucleic acid type, as a string.
	 * @param title   The title of the window.
	 */
	protected ModuleWindow( String acid, String title ) {

		// Initialize whether this window handles RNA or DNA.
		isRNA = ( acid.equals( "RNA" ) ) ? true : false;

		// Initialize the layout and constraints.
		layout = new GridBagLayout();
		setLayout( layout );
		constraints = new GridBagConstraints();

		// Create the progress bar dialog (which may or may not be used).
		progressBar = new JDialog();

		// Build the window.
		setTitles( title );
		makeInputControls();
	}

	/**
	 * Combine two file names into a new file name.
	 *
	 * @param file1   The first file name.
	 * @param file2   The second file name.
	 * @param ext     The extension of the new file name.
	 * @return        The combined file name.
	 */
	protected final String combineFileNames(
		String file1, String file2, String ext ) {
		String separator = File.separator;
		String name1 = new File( file1 ).getAbsolutePath();
		String name2 = new File( file2 ).getAbsolutePath();
		String dir = name1.substring( 0, name1.lastIndexOf( separator ) + 1 );
		name1 = name1.substring(
			name1.lastIndexOf( separator ) + 1,
			name1.lastIndexOf( "." ) );
		name2 = name2.substring(
			name2.lastIndexOf( separator ) + 1,
			name2.lastIndexOf( "." ) );
		return dir + name1 + "_" + name2 + "." + ext;
	}

	@Override
	protected void doActions( String command ) {

		// If the command comes from the "START" button, run the main
		// calculation in a new background thread.
		if( command.equals( "START" ) ) {
			new Thread() {
				public void run() { runMainCalculation(); }
			}.start();
		}

		// Otherwise, if the command comes from somewhere different, do one of
		// the window's specific actions.
		else { executeModuleAction( command ); }
	}

	/**
	 * Execute an action within a particular module.
	 *
	 * @param command   The command that signifies a particular action.
	 */
	protected abstract void executeModuleAction( String command );

	/**
	 * Handle finishing of a module.
	 * <br><br>
	 * Finishing of a module involves closing of the window and showing any
	 * error that may have occurred.
	 *
	 * @param text   The result of the calculation spawned by the module.
	 */
	protected final void finishModule( String text ) {
		this.dispose();
		progressBar.dispose();
		if( !text.equals( "" ) ) { dialogHandler.makeErrorDialog( text ); }
	}

	/**
	 * Get a particular input control from the content pane.
	 *
	 * @param index   The index of the control to get, one-indexed.
	 * @return        The input control.
	 */
	protected final JComponent getInputControl( int index ) {
		return getInputControl( getContentPane(), index );
	}

	/**
	 * Get a particular input control from a given container.
	 *
	 * @param hold    The container that holds the input control.
	 * @param index   The index of the control to get, one-indexed.
	 * @return        The input control.
	 */
	protected final JComponent getInputControl( Container hold, int index ) {
		return (JComponent)hold.getComponent( index - 1 );
	}

	/**
	 * Make any input controls necessary to get data from this window for a
	 * calculation.
	 */
	protected abstract void makeInputControls();

	/**
	 * Create the "START" button.
	 * <br><br>
	 * Note that this method simply makes the button, adds the action
	 * listener, applies standard padding, then places the button at the
	 * specified (X,Y) layout coordinate.
	 * <br><br>
	 * If more customization of the button itself or its positioning is
	 * needed, it should be done in subclasses.
	 *
	 * @param x   The X index at which to place the button in the layout.
	 * @param y   The Y index at which to place the button in the layout.
	 */
	protected final void makeStartButton( int x, int y ) {
		JButton startButton = new JButton( "START" );
		startButton.addActionListener( this );
		setPad( 15, 15 );
		setInsets( 10, 10, 10, 10 );
		setFillCenter();
		placeComponent( x, y, startButton );
	}

	/**
	 * Check to see if the module was initialized correctly.
	 * <br><br>
	 * Check what the result string retrieved from the module was, and return
	 * true only if it's the empty string. Otherwise, show the result string
	 * as an error and return false.
	 *
	 * @param result   The result of the data structure construction.
	 * @return         True if the data structure was initialized correctly,
	 *                 false if not.
	 */
	protected final boolean moduleInitialized( String result ) {
		if( result.equals( "" ) ) { return true; }
		else {
			dialogHandler.makeErrorDialog( result );
			return false;
		}
	}

	/**.
	 * Set a component in its place.
	 *
	 * @param x           The X coordinate of the grid.
	 * @param y           The Y coordinate of the grid.
	 * @param component   The component to set in place.
	 */
	protected final void placeComponent( int x, int y, JComponent component ) {
		constraints.gridx = x;
		constraints.gridy = y;
		add( component, constraints );
	}

	/**
	 * Replace the extension on a file name.
	 *
	 * @param file        The file name whose extension to replace.
	 * @param extension   The new extension.
	 * @return            The new file name with extension replaced.
	 */
	protected final String replaceExtension( String file, String extension ) {
		String root = file.substring( 0, file.lastIndexOf( "." ) );
		String newFile = root + "." + extension;
		return newFile;
	}

	/**
	 * Run the module calculation in the back end.
	 * <br><br>
	 * Data must be pulled from the module window before the back end
	 * calculation is run, so this data capture should be done in this method.
	 */
	protected abstract void runMainCalculation();

	/**
	 * Set the anchoring area for the GridBagConstraints to CENTER.
	 */
	protected final void setAnchorCenter() {
		constraints.anchor = GridBagConstraints.CENTER;
	}

	/**
	 * Set the anchoring area for the GridBagConstraints to NORTH.
	 */
	protected final void setAnchorNorth() {
		constraints.anchor = GridBagConstraints.NORTH;
	}

	/**
	 * Set the fill for the GridBagConstraints as CENTER.
	 */
	protected final void setFillCenter() {
		constraints.fill = GridBagConstraints.CENTER;
	}

	/**
	 * Set the fill for the GridBagConstraints as HORIZONTAL.
	 */
	protected final void setFillHorizontal() {
		constraints.fill = GridBagConstraints.HORIZONTAL;
	}

	/**
	 * Set the amount of space a component takes up in the grid.
	 *
	 * @param width    The width of the component.
	 * @param height   The height of the component.
	 */
	protected final void setGrid( int width, int height ) {
		constraints.gridwidth = width;
		constraints.gridheight = height;
	}

	/**
	 * Set the external padding of the constraints.
	 *
	 * @param top      The top padding.
	 * @param left     The left padding.
	 * @param bottom   The bottom padding.
	 * @param right    The right padding.
	 */
	protected final void setInsets(
		int top, int left, int bottom, int right ) {
		constraints.insets = new Insets( top, left, bottom, right );
	}

	/**
	 * Set the internal padding of the constraints.
	 *
	 * @param xPad   The padding in the X direction.
	 * @param yPad   The padding in the Y direction.
	 */
	protected final void setPad( int xPad, int yPad ) {
		constraints.ipadx = xPad;
		constraints.ipady = yPad;
	}

	/**
	 * Create a progress bar.
	 *
	 * @param determinate   True if the progress bar is determinate,
	 *                      false if not.
	 */
	private void showProgressBar( boolean determinate ) {

		// Build the progress bar.
		final JProgressBar bar = new JProgressBar();
		bar.setPreferredSize( new Dimension( 300, 50 ) );
		bar.setValue( 0 );

		// Add the progress bar to a panel.
		JPanel panel = new JPanel();
		panel.add( bar );
		new BorderBuilder().makeEqualBorder( 10, panel );

		// Add the progress bar panel to a dialog.
		progressBar.setTitle( "Calculation in Progress..." );
		progressBar.add( panel );
		progressBar.pack();
		progressBar.setVisible( true );

		// If the progress bar is indeterminate, explicitly set it as such.
		// This creates a simple entirely GUI-based dialog that just shows
		// that a calculation is happening. Otherwise, build a connection to
		// the back end to monitor and show explicit progress.
		if( !determinate ) { bar.setIndeterminate( true ); }
		else {

			// Build the SwingWorker that handles the background monitoring of
			// the progress bar. The SwingWorker runs a loop that, every half
			// second, checks the progress of the back end calculation and
			// updates the progress bar. The "done" method is called once the
			// calculation is done to dispose of the progress bar.
			SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
				public Void doInBackground() {
					try {
						setProgress( 0 );
						int progress = 0;

						while( progress < 100 ) {
							try { Thread.sleep( 500 ); }
							catch( InterruptedException ex ) {}

							progress = backend.getProgressNumber();
							setProgress( Math.min( progress, 100 ) );
						}
					} catch( Exception e ) { e.printStackTrace(); }

					return null;
				}

				public void done() { progressBar.dispose(); }
			};

			// Add a property change listener to the SwingWorker so the value
			// of the progress bar updates properly.
			worker.addPropertyChangeListener( new PropertyChangeListener() {
				public void propertyChange( PropertyChangeEvent e ) {
					if( e.getPropertyName() == "progress" ) {
						int progress = (Integer)e.getNewValue();
						bar.setValue( progress );
					}
				}
			});

			// Execute the progress bar through the SwingWorker.
			worker.execute();
		}
	}

	/**
	 * Create a determinate progress bar.
	 */
	protected void showProgressBarDeterminate() { showProgressBar( true ); }

	/**
	 * Create an indeterminate progress bar.
	 */
	protected void showProgressBarIndeterminate() { showProgressBar( false ); }
}
