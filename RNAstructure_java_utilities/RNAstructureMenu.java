/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

/**
 * A class that creates a menu.
 *
 * @author Jessica S. Reuter
 */
public abstract class RNAstructureMenu
	extends JMenu implements ActionListener {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param title   The title of the menu.
	 */
	protected RNAstructureMenu( String title ) {
		setText( title );
	}

	@Override
	public void actionPerformed( ActionEvent e ) {
		doMenuActions( e.getActionCommand() );
	}

	/**
	 * Add a check box menu item with rollover text.
	 *
	 * @param text       The menu item text.
	 * @param rollover   The rollover text.
	 */
	protected JCheckBoxMenuItem addCheckBoxMenuItem(
		String text, String rollover ) {
		JCheckBoxMenuItem item = new JCheckBoxMenuItem( text );
		item.setName( text );
		item.setSelected( false );
		item.setActionCommand( text );
		item.addActionListener( this );
		item.setToolTipText( rollover );
		add( item );
		return item;
	}

	/**
	 * Add a plain text menu item with rollover text.
	 *
	 * @param text       The menu item text.
	 * @param rollover   The rollover text.
	 */
	protected JMenuItem addMenuItem( String text, String rollover ) {
		JMenuItem item = new JMenuItem( text );
		item.setName( text );
		item.setActionCommand( text );
		item.addActionListener( this );
		item.setToolTipText( rollover );
		add( item );
		return item;
	}

	/**
	 * Add a plain text menu item with rollover text and a key stroke which
	 * can activate it.
	 *
	 * @param text       The menu item text.
	 * @param rollover   The rollover text.
	 * @param key        The key stroke.
	 */
	protected JMenuItem addMenuItem(
		String text, String rollover, char key ) {
		JMenuItem item = addMenuItem( text, rollover );
		int mask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();
		KeyStroke stroke = KeyStroke.getKeyStroke( key, mask );
		item.setAccelerator( stroke );
		return item;
	}

	/**
	 * Add a plain text menu item with rollover text and a key stroke which
	 * can activate it.
	 *
	 * @param text       The menu item text.
	 * @param rollover   The rollover text.
	 * @param key        The key stroke.
	 */
	protected JMenuItem addMenuItem(
		String text, String rollover, String key ) {
		JMenuItem item = addMenuItem( text, rollover );
		KeyStroke stroke = KeyStroke.getKeyStroke( key );
		item.setAccelerator( stroke );
		return item;
	}

	/**
	 * Do any actions specific to this menu.
	 *
	 * @param command   The command that signifies a particular action.
	 */
	protected abstract void doMenuActions( String command );

	/**
	 * Insert a plain text menu item with rollover text.
	 * <br><br>
	 * This method is used to add window-specific menu items and actions to an
	 * already existing menu. Note that any menu item added this way must have
	 * its actions defined in the window it acts upon.
	 *
	 * @param text       The menu item text.
	 * @param rollover   The rollover text.
	 * @param index      The index at which to insert the item, one-indexed.
	 * @param listener   The action listener attached to this item.
	 */
	public JMenuItem insertMenuItem(
		String text, String rollover, ActionListener listener, int index ) {
		JMenuItem item = new JMenuItem( text );
		item.setName( text );
		item.setActionCommand( text );
		item.addActionListener( listener );
		item.setToolTipText( rollover );
		insert( item, index - 1 );
		return item;
	}
}
