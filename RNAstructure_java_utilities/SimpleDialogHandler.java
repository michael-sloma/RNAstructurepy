/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.awt.BorderLayout;
import java.io.Serializable;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 * A class that handles simple dialogs used in the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 */
public class SimpleDialogHandler
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Show an error dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public void makeErrorDialog( String message ) {
		JOptionPane.showMessageDialog(
			null, "<html><center>" + message, "RNAstructure Error",
			JOptionPane.ERROR_MESSAGE );
	}

	/**
	 * Show an input dialog that gets a double within a certain range.
	 *
	 * @param field   The double input field.
	 * @return        The double selected from this field.
	 */
	public Double makeInputDialog( NumberField.DoubleField field ) {
		JPanel panel = new JPanel( new BorderLayout() );
		panel.add( new JLabel( field.getName() ), BorderLayout.NORTH );
		panel.add( field );
		JOptionPane.showOptionDialog(
			null, panel, "RNAstructure Input", JOptionPane.OK_CANCEL_OPTION,
			JOptionPane.INFORMATION_MESSAGE, null, null, null );
		return field.getValue();
	}

	/**
	 * Show an input dialog that gets an integer within a certain range.
	 *
	 * @param field   The integer input field.
	 * @return        The integer selected from this field.
	 */
	public Integer makeInputDialog( NumberField.IntegerField field ) {
		JPanel panel = new JPanel( new BorderLayout() );
		panel.add( new JLabel( field.getName() ), BorderLayout.NORTH );
		panel.add( field );
		JOptionPane.showOptionDialog(
			null, panel, "RNAstructure Input", JOptionPane.OK_CANCEL_OPTION,
			JOptionPane.INFORMATION_MESSAGE, null, null, null );
		return field.getValue();
	}

	/**
	 * Show a message dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public void makeMessageDialog( Object message ) {
		JOptionPane.showMessageDialog(
			null, message, "RNAstructure", JOptionPane.PLAIN_MESSAGE );
	}

	/**
	 * Show a message dialog.
	 *
	 * @param msg   The message on the dialog.
	 */
	public void makeMessageDialog( String msg ) {
		StringBuilder builder = new StringBuilder( "<html><center>" + msg );
		makeMessageDialog( builder );
	}

	/**
	 * Show a dialog that gives the user one of three choices to make.
	 *
	 * @param message   The message on the dialog.
	 * @return          The response: "YES", "NO", or "CANCEL".
	 */
	public String makeThreeChoicesDialog( Object message ) {
		int choice = JOptionPane.showConfirmDialog(
			null, message, "RNAstructure",
			JOptionPane.YES_NO_CANCEL_OPTION );
		if( choice == JOptionPane.YES_OPTION ) { return "YES"; }
		else if( choice == JOptionPane.NO_OPTION ) { return "NO"; }
		else { return "CANCEL"; }
	}

	/**
	 * Show a dialog that gives the user one of three choices to make.
	 *
	 * @param msg   The message on the dialog.
	 * @return      The response: "YES", "NO", or "CANCEL".
	 */
	public String makeThreeChoicesDialog( String msg ) {
		StringBuilder builder = new StringBuilder( "<html><center>" + msg );
		return makeThreeChoicesDialog( builder );
	}

	/**
	 * Show a dialog that gives the user one of two choices to make.
	 *
	 * @param message   The message on the dialog.
	 * @return          The response: "OK" or "CANCEL".
	 */
	public String makeTwoChoicesDialog( Object message ) {
		int choice = JOptionPane.showConfirmDialog(
			null, message, "RNAstructure", JOptionPane.OK_CANCEL_OPTION );
		if( choice == JOptionPane.OK_OPTION ) { return "OK"; }
		else { return "CANCEL"; }
	}

	/**
	 * Show a dialog that gives the user one of two choices to make.
	 *
	 * @param msg   The message on the dialog.
	 * @return      The response: "OK" or "CANCEL".
	 */
	public String makeTwoChoicesDialog( String msg ) {
		StringBuilder builder = new StringBuilder( "<html><center>" + msg );
		return makeTwoChoicesDialog( builder );
	}

	/**
	 * Show a warning dialog.
	 *
	 * @param message   The message on the dialog.
	 */
	public void makeWarningDialog( String message ) {
		JOptionPane.showMessageDialog(
			null, "<html><center>" + message, "RNAstructure",
			JOptionPane.WARNING_MESSAGE );
	}
}
