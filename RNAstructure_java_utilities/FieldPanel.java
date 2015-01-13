/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * A class that creates a panel which holds aligned input fields.
 *
 * @author Jessica S. Reuter
 */
public class FieldPanel
	extends JPanel {
	private static final long serialVersionUID = 20120802;

	/**
	 * The components on this panel (buttons or labels).
	 */
	protected JComponent[] components;

	/**
	 * The specified width of the panel.
	 */
	private int specifiedPanelWidth;

	/**
	 * The array of text fields.
	 */
	protected JTextField[] textFields;

	/**
	 * Constructor.
	 *
	 * @param fields   The array of text fields in this panel.
	 */
	protected FieldPanel( JTextField[] fields ) {
		setLayout( new GridLayout( 0, 1 ) );
		textFields = fields;
		components = new JComponent[textFields.length];
		specifiedPanelWidth = 0;
	}

	/**
	 * Build the panel.
	 */
	public void makePanel() {

		// Determine the border that should go around the panel.
		int border = 10;

		// Determine the width of the widest component associated with the
		// fields, so all components can be that width and align.
		int maxWidth = 0;
		int numFields = textFields.length;
		for( int i = 1; i <= numFields; i++ ) {
			int nextWidth = components[i-1].getPreferredSize().width;
			maxWidth = Math.max( maxWidth, nextWidth );
		}

		// For each component in the array, place it next to its field.
		for( int i = 1; i <= numFields; i++ ) {

			// Set the preferred width of the component associated with the
			// data field.
			JComponent component = components[i-1];
			int height = component.getPreferredSize().height;
			component.setPreferredSize( new Dimension( maxWidth, height ) );

			// Create the panel that holds the new component and text field,
			// then add it to the input panel.
			JPanel panel = new JPanel( new BorderLayout() );
			panel.add( component, BorderLayout.WEST );
			panel.add( textFields[i-1] );
			add( panel );
		}

		// Set the full panel size.
		int preferredHeight =
			( new JTextField().getPreferredSize().height * numFields ) +
			( 2 * border ); 
		Dimension preferredDimension =
			new Dimension( specifiedPanelWidth, preferredHeight );
		setPreferredSize( preferredDimension );
		setMinimumSize( preferredDimension );
		setMaximumSize( preferredDimension );
		setSize( preferredDimension );

		// Set a border around the panel.
		new BorderBuilder().makeEqualBorder( border, this );
	}

	/**
	 * Set the width of the panel.
	 */
	public void setPanelWidth( int width ) { specifiedPanelWidth = width; }

	/**
	 * An inner class that handles a panel which contains file fields.
	 *
	 * @author Jessica S. Reuter
	 */
	public static class FilePanel
		extends FieldPanel {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 * <br><br>
		 * Create a panel of file fields attached to buttons.
		 *
		 * @param listener   The action listener attached to the buttons.
		 * @param fields     The file fields.
		 */
		public FilePanel(
			final ActionListener listener, final FileField... fields ) {
			super( fields );

			// For each text field, create an element on the panel.
			for( int i = 1; i <= fields.length; i++ ) {
				final int index = i;

				// Create the file button.
				JButton button = new JButton( fields[i-1].getName() );
				components[i-1] = button;

				// Add an action listener to the button.
				button.addActionListener( new ActionListener() {
					public void actionPerformed( ActionEvent e ) {

						// Set the button's ability to select files as true.
						boolean canSelect = true;

						// If the button is not the first, check to make sure
						// all previous fields have been filled.
						// If not, set its ability to select as false.
						if( index != 1 ) {
							for( int i = 1; i <= index - 1; i++ ) {
								if( getFile( i ).equals( "" ) ) {
									String message =
										fields[index-1].getName() +
										" cannot be selected at this time.";
									SimpleDialogHandler handler =
										new SimpleDialogHandler();
									handler.makeErrorDialog( message );
									canSelect = false;
									break;
								}
							}
						}

						// If the button has the ability to select files, do
						// the file selection based on the listener passed in
						// to the function.
						if( canSelect ) { listener.actionPerformed( e ); }
					}
				});
			}
		}

		/**
		 * Get a file name from a particular input field.
		 *
		 * @param index   The text field whose contents to get, one-indexed.
		 * @return        The file name.
		 */
		public String getFile( int index ) {
			return textFields[index-1].getText().trim();
		}

		/**
		 * Get whether any fields are empty; an empty field is an error.
		 * 
		 * @return   Whether an error happened.
		 */
		public boolean isError() {
			for( JTextField element: textFields ) {
				if( element.getText().trim().equals( "" ) ) {
					String msg = element.getName() + " must not be empty.";
					new SimpleDialogHandler().makeErrorDialog( msg );
					return true;
				}
			}
			return false;
		}

		/**
		 * Set a file name in a particular input field.
		 *
		 * @param index   The text field whose contents to set, one-indexed.
		 * @param file    The file name to set in the field.
		 */
		public void setFile( int index, String file ) {
			textFields[index-1].setText( file );
		}
	}

	/**
	 * An inner class that handles a panel which contains numeric data fields.
	 *
	 * @author Jessica S. Reuter
	 */
	public static class NumberPanel
		extends FieldPanel {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 * <br><br>
		 * Create a panel of numeric data fields attached to labels.
		 *
		 * @param fields   The data fields.
		 */
		public NumberPanel( NumberField... fields ) {
			super( fields );
			for( int i = 1; i <= fields.length; i++ ) {
				JLabel label = new JLabel( fields[i-1].getName() );
				label.setHorizontalAlignment( JLabel.LEFT );
				new BorderBuilder().makeRightBorder( 2, label );
				components[i-1] = label;
			}
		}

		/**
		 * Get a  particular input field.
		 *
		 * @param index   The text field to get, one-indexed.
		 * @return        The input field.
		 */
		public JTextField getField( int index ) {
			return textFields[index-1];
		}
	}
}
