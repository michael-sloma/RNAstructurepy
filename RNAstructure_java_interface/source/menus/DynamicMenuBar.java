/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.menus;

import javax.swing.JMenuBar;

import RNAstructure_java_interface.source.windows.InternalWindow;

/**
 * A class that creates a menu bar.
 *
 * @author Jessica S. Reuter
 */
public class DynamicMenuBar
	extends JMenuBar {
	private static final long serialVersionUID = 20120802;

	/**
	 * Constructor.
	 *
	 * @param menus   The variable menus this bar contains.
	 * @param win     The window this menu bar is connected to.
	 */
	public DynamicMenuBar( InternalWindow win, RolloverMenu... menus ) {

		// Create the File menu.
		FileMenu file = new FileMenu( win );
		file.convertToolTips();
		add( file );

		// Create the RNA menu.
		NucleicAcidMenu rna = NucleicAcidMenu.createRNA();
		rna.convertToolTips();
		add( rna );

		// Create the DNA menu.
		NucleicAcidMenu dna = NucleicAcidMenu.createDNA();
		dna.convertToolTips();
		add( dna );

		// Create the variable menus.
		if( menus != null ) {
			for( RolloverMenu menu: menus ) {
				menu.setEnabled( false );
				menu.convertToolTips();
				add( menu );
			}
		}

		// Create the Help menu.
		HelpMenu help = new HelpMenu();
		help.convertToolTips();
		add( help );
	}

	/**
	 * Enable the menus on this menu bar.
	 */
	public void enableMenus() {
		int menus = getMenuCount();
		for( int i = 1; i <= menus; i++ ) {
			getMenu( i - 1 ).setEnabled( true );
		}
	}
}
