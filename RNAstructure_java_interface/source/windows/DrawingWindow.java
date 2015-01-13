/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.swing.JMenuBar;

import RNAstructure_java_drawing.source.dialogs.ImageDialog;
import RNAstructure_java_drawing.source.dialogs.PlotDialog;
import RNAstructure_java_drawing.source.dialogs.StructureDialog;
import RNAstructure_java_interface.source.menus.DynamicMenuBar;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.MacChecker;
import RNAstructure_java_utilities.RNAstructureMenu;
import RNAstructure_java_utilities.ScrollerPane;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for displaying a drawing window.
 * <br><br>
 * This class is a wrapper for the self-contained drawing framework dialogs,
 * which essentially converts them into internal frames but does not take over
 * or duplicate any of their functions.
 * <br><br>
 * Most drawings, menus, or listeners are transplanted from the drawing
 * framework dialogs.
 *
 * @author Jessica S. Reuter
 */
public class DrawingWindow
	extends InternalWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The image dialog connected to this window.
	 */
	private ImageDialog dialog;

	/**
	 * True if an error occurred building this window, false if not.
	 */
	private boolean error;

	/**
	 * The strand being drawn.
	 */
	private int strand;

	/**
	 * File Constructor.
	 *
	 * @param file   The file whose data should be drawn on the window.
	 */
	public DrawingWindow( String file ) {
		this( file, 1 );
	}

	/**
	 * Strand Constructor.
	 * 
	 * @param file     The dot plot data file to draw.
	 * @param strand   The strand from the data file to draw.
	 */
	public DrawingWindow( String file, int strand ) {

		// Set the title of this window and the main frame to reflect the
		// creation of a new drawing. Also, set the frame to be resizable.
		setTitles( file );
		setResizable( true );

		// Set the strand.
		this.strand = strand;

		// If the file is a CT file, get information from a structure dialog.
		// Otherwise, get information from a plot dialog.
		if( file.endsWith( "ct" ) ) { dialog = new StructureDialog( file ); }
		else { dialog = new PlotDialog( file, strand ); }
		getImageData();

		// Add a component listener to the legend so every time it changes, the
		// window is revalidated.
		ComponentListener cListener = new ComponentListener() {

			@Override
			public void componentHidden( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentMoved( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentResized( ComponentEvent e ) { doRevalidate(); }

			@Override
			public void componentShown( ComponentEvent e ) { doRevalidate(); }
		};
		getContentPane().getComponent( 2 ).addComponentListener( cListener );

		// If no error occurred and the dialog is a structure dialog, grab its
		// structure movement key listener.
		if( ( error == false ) && ( file.endsWith( "ct" ) ) ) {
			int index = dialog.getKeyListeners().length - 1;
			addKeyListener( dialog.getKeyListeners()[index] );
		}

		// If no error occurred and the dialog is a probability plot dialog,
		// replace the listener on the "Write Probable Structures" menu item
		// with this class, so structures can pop up in DrawingWindows instead
		// of StructureDialogs.
		if( ( error == false ) && ( file.endsWith( "pfs" ) ) ) {
			ActionListener aListener =
				menuBar.getMenu( 4 ).getItem( 3 ).getActionListeners()[0];
			menuBar.getMenu( 4 ).getItem( 3 ).removeActionListener( aListener );
			menuBar.getMenu( 4 ).getItem( 3 ).addActionListener( this );
		}
	}

	/**
	 * {@inheritDoc}
	 * <br><br>
	 * Note that most of the actions for this window are handled by the menus
	 * in the self-contained drawing framework.
	 */
	@Override
	protected void doActions( String command ) {

		// If the command is to write a probable structures file, do so.
		if( command.equals( "Write Probable Structures File" ) ) {

			// Choose and write a structure file.
			String file = FileChooser.doSave().getCT();
			boolean drawOK = false;
			if( !file.equals( "" ) ) {
				drawOK = ((PlotDialog)dialog).writeStructuresFile(
					file, false );
			}

			// If the user wants to draw structures, draw them.
			if( drawOK ) {
				String query = "Do you want to draw structures?";
				String doDraw = new SimpleDialogHandler()
					.makeTwoChoicesDialog( query );
				if( doDraw.equals( "OK" ) ) {
					DrawingWindow window = new DrawingWindow( file );
					if( window.isError() == false ) { window.viewWindow(); }
				}
			}
		}
	}

	/**                                                                                                                                          
	 * Revalidate the drawing window.                                                                                                            
	 */
	private void doRevalidate() {

		// Get the content pane.
		Container pane = getContentPane();

		// Get the main window components.
		Container label = (Container)pane.getComponent( 0 );
		ScrollerPane mainScroll = (ScrollerPane)pane.getComponent( 1 );
		ScrollerPane legendScroll = (ScrollerPane)pane.getComponent( 2 );
		Container legend = (Container)(legendScroll.getViewport().getView());

		// Calculate the container pane size and set it.
		Dimension size = mainScroll.getSize();
		size.height += label.getPreferredSize().height;
		if( legend.getComponents().length != 0 ) {
			size.height += legendScroll.getPreferredSize().height;
		}
		getContentPane().setPreferredSize( size );

		// Repaint the window.                                                                                                                  
		repaint();
	}

	/**
	 * Get the data from a dialog to put in this window.
	 */
	private void getImageData() {

		// If the dialog encountered an error while being built,
		// set the error flag and return.
		if( dialog.isError() ) {
			error = true;
			return;
		}

		// Set the dialog's content pane as the frame's content pane.
		setContentPane( dialog.getContentPane() );

		// Convert the dialog's menus into rollover menus and set them in the
		// dynamic menu bar.
		JMenuBar bar = dialog.getJMenuBar();
		RolloverMenu[] convertedMenus = new RolloverMenu[bar.getMenuCount()];
		for( int i = 1; i <= convertedMenus.length; i++ ) {
			RNAstructureMenu menu = (RNAstructureMenu)bar.getMenu( i - 1 );
			convertedMenus[i-1] = new RolloverMenu( menu );
		}
		menuBar = new DynamicMenuBar( this, convertedMenus );
		menuBar.enableMenus();

		// Grab the dialog's zooming key listener.
		addKeyListener( dialog.getKeyListeners()[0] );

		// Add the dialog's scroll bar listener, if the OS is a type of Mac.
		if( MacChecker.isMac() ) {
			addKeyListener( dialog.getKeyListeners()[1] );
		}
	}

	/**
	 * Get whether this window encountered an error.
	 *
	 * @return   True if an error occurred, false if not.
	 */
	public boolean isError() { return error; }

	/**
	 * {@inheritDoc}
	 * <br><br>
	 * Note that this version of the function does nothing, since this window
	 * is only a wrapper for a self-contained drawing dialog.
	 * <br><br>
	 * Menus are added dynamically in the drawing window class, because this
	 * class can hold multiple types of drawings with different menus.
	 */
	@Override
	protected RolloverMenu[] setMenus() { return null; }
}
