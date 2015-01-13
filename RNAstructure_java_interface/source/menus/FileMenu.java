/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.menus;

import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JButton;
import RNAstructure_java_interface.source.startup.ApplicationRootFrame;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_interface.source.windows.DrawingWindow;
import RNAstructure_java_interface.source.windows.DynalignRefoldWindow;
import RNAstructure_java_interface.source.windows.InternalWindow;
import RNAstructure_java_interface.source.windows.OligoScreenWindow;
import RNAstructure_java_interface.source.windows.RefoldWindow;
import RNAstructure_java_interface.source.windows.SequenceDisplayWindow;

/**
 * A class that creates a "File" menu.
 * <br><br>
 * Most File menus are identical, but items can be added to or subtracted from
 * them in specific contexts.
 *
 * @author Jessica S. Reuter
 */
public class FileMenu
	extends RolloverMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param window   The window this menu is attached to.
	 */
	public FileMenu( InternalWindow window ) {
		super( "File" );
		addMenuItem( "New Sequence", "Create a new sequence.", 'N' );
		addMenuItem( "Open Sequence", "Open an existing sequence.", 'O' );
		addSeparator();
		addMenuItem( "OligoScreen",
			"Calculate thermodynamic parameters for a set of " +
			"oligonucleotides." );
		addSeparator();
		addMenuItem( "Draw",
			"Draw a secondary structure." );
		addMenuItem( "Dot Plot",
			"Display the energy dot plot for a sequence that was previously " +
			"folded." );
		addMenuItem( "Dot Plot Partition Function",
			"Display base pairing probabilities for a previously calculated " +
			"sequence." );
		addMenuItem( "Dot Plot Dynalign",
			"Generate a Dynalign dot plot for two sequences." );
		addMenuItem( "Dot Plot From Text File",
			"Draw a dot plot from a text file." );
		addSeparator();
		addMenuItem( "Refold From Save File",
			"Refold a sequence from its save file." );
		addMenuItem( "Refold From Dynalign Save File",
			"Refold from a Dynalign calculation." );
		addSeparator();
		addMenuItem( "Exit", "Exit the RNAstructure application." );
	}

	@Override
	protected void doMenuActions( String command ) {

		// If the command is a drawing command, prepare a drawing window.
		if( command.startsWith( "D" ) ) {

			// Show an energy dot plot window, if possible.
			if( command.equals( "Dot Plot" ) ) {
				String file = FileChooser.doOpen().getFolding();
				if( !file.equals( "" ) ) {
					DrawingWindow window = new DrawingWindow( file );
					if( window.isError() == false ) { window.viewWindow(); }
				}
			}

			// Show a Dynalign dot plot window, if possible.
			else if( command.equals( "Dot Plot Dynalign" ) ) {
				String file = FileChooser.doOpen().getDynalign();
				if( !file.equals( "" ) ) {
					DrawingWindow window1 = new DrawingWindow( file, 1 );
					DrawingWindow window2 = new DrawingWindow( file, 2 );
					boolean validWindows =
						( window1.isError() == false ) &&
						( window2.isError() == false );
					if( validWindows ) {
						window1.viewWindow();
						window2.viewWindow();
					}
				}
			}

			// Show a text dot plot window, if possible.
			else if( command.equals( "Dot Plot From Text File" ) ) {
				String file = FileChooser.doOpen().getDotPlot();
				if( !file.equals( "" ) ) {
					DrawingWindow window = new DrawingWindow( file );
					if( window.isError() == false ) { window.viewWindow(); }
				}
			}

			// Show a probability dot plot window, if possible.
			else if( command.equals( "Dot Plot Partition Function" ) ) {
				String file = FileChooser.doOpen().getPartition();
				if( !file.equals( "" ) ) {
					DrawingWindow window = new DrawingWindow( file );
					if( window.isError() == false ) { window.viewWindow(); }
				}
			}

			// Initialize a structure window, if possible.
			else if( command.equals( "Draw" ) ) {
				String file = FileChooser.doOpen().getCT();
				if( !file.equals( "" ) ) {
					DrawingWindow window = new DrawingWindow( file );
					if( window.isError() == false ) { window.viewWindow(); }
				}
			}
		}

		// If the command is to exit the application, do so.
		else if( command.equals( "Exit" ) ) { System.exit( 0 ); }

		// If the command is to build a new sequence, show a blank sequence
		// display window.
		else if( command.equals( "New Sequence" ) ) {
			new SequenceDisplayWindow().viewWindow();
		}

		// If the command is to run OligoScreen, show an OligoScreen window.
		else if( command.equals( "OligoScreen" ) ) {
			new OligoScreenWindow().viewWindow();
		}

		// If the command is to open an existing sequence, show a filled
		// sequence display window.
		else if( command.equals( "Open Sequence" ) ) {
			String file = FileChooser.doOpen().getSequenceExtended();
			if( !file.equals( "" ) ) {
				new SequenceDisplayWindow( file ).viewWindow();
			}
		}

		// If the command is to refold from a Dynalign save file, open a
		// Dynalign refolding window.
		else if( command.equals( "Refold From Dynalign Save File" ) ) {
			new DynalignRefoldWindow().viewWindow();
		}

		// If the command is to refold from a folding save file, open a
		// refolding window.
		else if( command.equals( "Refold From Save File" ) ) {
			new RefoldWindow().viewWindow();
		}
	}
}
