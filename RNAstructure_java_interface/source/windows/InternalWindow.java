/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import javax.swing.ImageIcon;
import javax.swing.JDesktopPane;
import javax.swing.JInternalFrame;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;

import RNAstructure_java_interface.source.menus.DynamicMenuBar;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_interface.source.startup.ApplicationRootFrame;
import RNAstructure_java_interface.source.startup.RNAstructureBackendCalculator;
import RNAstructure_java_utilities.MacChecker;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class that holds internal frames for the RNAstructure GUI.
 * <br><br>
 * Every frame is connected to the back end through an individual proxy.
 *
 * @author Jessica S. Reuter
 */
public abstract class InternalWindow
	extends JInternalFrame implements ActionListener {
	private static final long serialVersionUID = 20120802;

	/**
	 * The back end calculator that's linked to this window.
	 */
	protected final RNAstructureBackendCalculator backend;

	/**
	 * A dialog handler to show any simple dialog messages that pop up from
	 * this window.
	 */
	protected final SimpleDialogHandler dialogHandler;

	/**
	 * The menu bar that's attached to this window.
	 */
	protected DynamicMenuBar menuBar;

	/**
	 * Constructor.
	 */
	protected InternalWindow() {

		// Initialize the back end calculator and the dialog handler.
		backend = new RNAstructureBackendCalculator();
		dialogHandler = new SimpleDialogHandler();

		// Set the frame defaults.
		setDefaultCloseOperation( JInternalFrame.DISPOSE_ON_CLOSE );
		setClosable( true );
		setFocusable( true );
		setIconifiable( false );
		setMaximizable( false );
		setResizable( false );

		// Add the window to the base frame.
		JDesktopPane desktop =
			(JDesktopPane)ApplicationRootFrame.getFrame()
			.getContentPane().getComponent( 1 );
		desktop.add( this );

		// Set the frame icon, if the look and feel allows for it.
		// If the OS is anything but a Mac, a frame icon can be set.
		if( MacChecker.notMac() ) {
		    ImageIcon icon =
				new ImageIcon( ApplicationRootFrame.getFrame().getIconImage() );
		    setFrameIcon( icon );
		}

		// Set the menu bar and menus for this window.
		menuBar = new DynamicMenuBar( this, setMenus() );

		// Add a focus listener so the window sets this window's associated
		// title and menu bar in the main frame.
		final JInternalFrame internal = this;
		addFocusListener( new FocusAdapter() {
			public void focusGained( FocusEvent e ) {
				ApplicationRootFrame frame =
					(ApplicationRootFrame)getDesktopPane()
					.getTopLevelAncestor();
				setTitles( getTitle() );
				frame.setJMenuBar( menuBar );
				frame.repaint();
				JDesktopPane desktop =
					(JDesktopPane)ApplicationRootFrame.getFrame()
					.getContentPane().getComponent( 1 );
				desktop.getDesktopManager().activateFrame( internal );
			}
		});

		// Add a window closing listener that resets the main frame to its
		// original state when no more windows are open.
		addInternalFrameListener( new InternalFrameAdapter() {
			public void internalFrameClosed( InternalFrameEvent e ) {
				ApplicationRootFrame frame = ApplicationRootFrame.getFrame();
				if( frame.getMostRecentFrame() == null ) {
					frame.resetFrame();
				}
			}
		});
	}

	@Override
	public void actionPerformed( ActionEvent e ) {
		doActions( e.getActionCommand() );
	}

	/**
	 * Do any actions specific to this window.
	 *
	 * @param command   The command that signifies a particular action.
	 */
	protected abstract void doActions( String command );

	/**
	 * Create the array of variable menus used by this window.
	 *
	 * @return   The array of variable menus created for this window.
	 */
	protected abstract RolloverMenu[] setMenus();

	/**
	 * Set the title of both this window and the main window.
	 *
	 * @param title   The new title.
	 */
	protected void setTitles( String title ) {
		setTitle( title );
		ApplicationRootFrame.getFrame().setTitle( "RNAstructure -- " + title );
	}

	/**
	 * View the window.
	 */
	public void viewWindow() {
		pack();
		setLocation( 0, 0 );
		setVisible( true );
		menuBar.revalidate();
		menuBar.repaint();
	}
}
