/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_interface.source.menus;

import java.awt.Component;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.Serializable;

import javax.swing.JMenuItem;

import RNAstructure_java_interface.source.startup.ApplicationRootFrame;
import RNAstructure_java_utilities.RNAstructureMenu;

/**
 * A class that creates a menu for the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 */
public class RolloverMenu
	extends RNAstructureMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * RNAstructureMenu Conversion Constructor.
	 *
	 * @param oldMenu   The menu to convert.
	 */
	public RolloverMenu( RNAstructureMenu oldMenu ) {
		super( oldMenu.getText() );
		Component[] components = oldMenu.getMenuComponents();
		for( Component component: components ) {
			if( component instanceof JMenuItem ) {
				add( (JMenuItem)component );
			}
		}
		convertToolTips();
	}

	/**
	 * Title Constructor.
	 *
	 * @param title   The title of the menu.
	 */
	protected RolloverMenu( String title ) {
		super( title );
	}

	/**
	 * Convert tooltips to rollover text.
	 */
	public void convertToolTips() {
		Component[] components = getMenuComponents();
		for( Component component: components ) {
			if( component instanceof JMenuItem ) {
				String tooltip = ((JMenuItem)component).getToolTipText();
				((JMenuItem)component).setToolTipText( null );
				((JMenuItem)component).addMouseListener(
					new RolloverListener( tooltip ) );
			}
		}
	}

	@Override
	protected void doMenuActions( String command ) { /* Do nothing. */ }

	/**
	 * An inner class that creates a rollover listener for menus in the
	 * RNAstructure GUI.
	 * <br><br>
	 * The rollover text is placed in the message bar at the bottom of the main
	 * application frame, whenever a user interacts with the menu item this
	 * listener is attached to, and removed to make way for the default message
	 * text when the user stops interacting with the menu item.
	 *
	 * @author Jessica S. Reuter
	 */
	public class RolloverListener
		implements MouseListener, Serializable {
		private static final long serialVersionUID = 20120802;

		/**
		 * The rollover text for this listener.
		 */
		private String rollover;

		/**
		 * Constructor.
		 *
		 * @param text   The rollover text.
		 */
		public RolloverListener( String text ) { rollover = text; }

		@Override
		public void mouseClicked( MouseEvent e ) { reset(); }

		@Override
		public void mouseEntered( MouseEvent e ) { set(); }

		@Override
		public void mouseExited( MouseEvent e ) { reset(); }

		@Override
		public void mousePressed( MouseEvent e ) { reset(); }

		@Override
		public void mouseReleased( MouseEvent e ) { reset(); }

		/**
		 * Reset the main frame message bar text to its default.
		 */
		private void reset() {
			ApplicationRootFrame.getFrame().resetInfoLabel();
		}

		/**
		 * Set the main frame message bar text to this listener's specific
		 * rollover text.
		 */
		private void set() {
			ApplicationRootFrame.getFrame().setInfoLabel( rollover );
		}
	}
}
