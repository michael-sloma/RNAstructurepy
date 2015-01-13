/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.menus;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JPanel;

import RNAstructure_java_interface.source.startup.RNAstructureBackendCalculator;
import RNAstructure_java_utilities.FieldPanel.FilePanel;
import RNAstructure_java_utilities.FieldPanel.NumberPanel;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.FileField;
import RNAstructure_java_utilities.NumberField;
import RNAstructure_java_utilities.NumberField.DoubleField;
import RNAstructure_java_utilities.NumberField.IntegerField;
import RNAstructure_java_utilities.RadioButtonPanel;
import RNAstructure_java_utilities.SimpleDialogHandler;
import RNAstructure_java_utilities.ValueSelectionDialog;

/**
 * A class that creates a folding constraints menu.
 *
 * @author Jessica S. Reuter
 */
public class ConstraintsMenu
	extends RolloverMenu {
	private static final long serialVersionUID = 20120802;

	/**
	 * The back end calculator that's linked to this menu.
	 */
	private RNAstructureBackendCalculator backend;

	/**
	 * The strand this constraints menu deals with.
	 */
	private int str = 0;

	/**
	 * General Constructor.
	 *
	 * @param calc   The back end calculator connected to this menu.
	 */
	public ConstraintsMenu( RNAstructureBackendCalculator calc ) {
		super( "Force" );
		backend = calc;
	}

	/**
	 * Add a menu section of general constraints commands.
	 */
	public void addGeneralSection() {
		addMenuItem( "Force Base Pairs",
			"Force two bases to be paired in a folded structure." );
		addMenuItem( "Prohibit Base Pairs",
			"Mandate that certain bases be unable to pair in a structure." );
		addMenuItem( "Chemical Modification",
			"Modify specified bases chemically prior to folding." );
		addMenuItem( "FMN Cleavage", "Cleave a strand of interest." );
		addMenuItem( "Single Stranded",
			"Force particular regions to be single stranded." );
		addMenuItem( "Double Stranded",
			"Force particular regions to be double stranded." );
		addSeparator();
	}

	/**
	 * Add a menu section that handles maximum pairing distance.
	 */
	public void addMaxPairingDistanceSection() {
		addMenuItem( "Maximum Pairing Distance",
			"Limit the maximum distance allowed between paired bases." );
		addSeparator();
	}

	/**
	 * Add a menu section that handles SHAPE constraints.
	 */
	public void addSHAPESection() {
		String[] types = { "Hard", "Pseudo-Energy" };
		for( String type: types ) {
			addMenuItem( "Read SHAPE Reactivity -- " + type + " Constraints",
				"Input constraints on SHAPE data to be used in folding." );
		}
		addSeparator();
	}

	/**
	 * Add a menu section that handles saving and restoring constraints.
	 */
	public void addSaveRestoreSection() {
		addMenuItem( "Save Constraints",
			"Save folding constraints to a save file." );
		addMenuItem( "Restore Constraints",
			"Get a set of folding constraints from a save file." );
		addSeparator();
	}

	/**
	 * Add a menu section that handles showing and resetting constraints.
	 */
	public void addShowResetSection() {
		addMenuItem( "Show Current Constraints",
			"View constraints that are currently applied to this sequence." );
		addMenuItem( "Reset Current Constraints",
			"Reset the constraints currently applied to this sequence." );
		addSeparator();
	}

	/**
	 * Build a complete menu that handles if unimolecular pairs are allowed.
	 */
	public void buildUnimolecularMenu() {
		addCheckBoxMenuItem( "Forbid Unimolecular Pairs",
			"Forbid pairs from forming between bases of the same molecule." );
	}

	/**
	 * Build a complete menu that handles Dynalign alignment constraints.
	 */
	public void buildDynalignAlignmentMenu() {
		setText( "Constraints for Alignment" );
		addMenuItem( "Force Alignment",
			"Force the alignment of particular bases." );
		addSeparator();
		addMenuItem( "Show Current Alignment Constraints",
			"Show constraints set for this particular alignment." );
		addMenuItem( "Reset Alignment Constraints",
			"Reset the constraints that are applied to this alignment." );
		addSeparator();
		addMenuItem( "Save Alignment", "Save the alignment to a save file." );
		addMenuItem( "Restore Alignment From File",
			"Get an alignment from a save file." );
	}

	/**
	 * Build a complete menu that handles maximum loop size.
	 */
	public void buildMaxLoopMenu() {
		setText( "Maximum Loop" );
		addMenuItem( "Set Maximum Loop Size",
			"Set the maximum number of unpaired nucleotides allowed." );
	}

	/**
	 * Build a complete menu that handles temperature.
	 *
	 * @return   The temperature menu.
	 */
	public ConstraintsMenu buildTemperatureMenu() {
		setText( "Temperature" );
		addMenuItem( "Set Temperature",
			"Set the temperature at which calculations occur." );
		return this;
	}

	@Override
	@SuppressWarnings("serial")
	protected void doMenuActions( String command ) {

		// Create a dialog handler to use for dialog generation.
		SimpleDialogHandler handler = new SimpleDialogHandler();

		// If the command handles alignment constraints, handle a particular
		// action referring to them, then return out of this method.
		if( command.contains( "Alignment" ) ) {

			// If the command is to force an alignment, show a dialog that
			// forces two nucleotides to align.
			if( command.startsWith( "Force" ) ) {
				String alignHeaders = "Base in Sequence 1;Base in Sequence 2";
				new NucOrHelixDialog( command, alignHeaders ) {
					public String setConstraint( String data ) {
						String[] values = data.split( " " );
						Integer nuc1 = Integer.parseInt( values[0] );
						Integer nuc2 = Integer.parseInt( values[1] );
						return backend.setDynalignAlignmentConstraint(
							nuc1, nuc2 );
					}
				};
			}

			// If the command is to reset alignment constraints, do so.
			else if( command.startsWith( "Reset" ) ) {
				boolean clear = showDeletionConfirmation();
				if( clear ) { backend.clearDynalignAlignmentConstraints(); }
			}

			// If the command is to read a constraints file, do so.
			else if( command.startsWith( "Restore" ) ) {
				String file = FileChooser.doOpen().getConstraints();
				String result =
					backend.readDynalignAlignmentConstraintsFile( file );
				if( !result.equals( "" ) ) {
					handler.makeErrorDialog( result );
				}
			}

			// If the command is to write a constraints file, do so.
			else if( command.startsWith( "Save" ) ) {
				String file = FileChooser.doSave().getConstraints();
				backend.writeDynalignAlignmentConstraintsFile( file );
				handler.makeMessageDialog(
					"Dynalign alignment constraints file written." );
			}

			// If the command is to show constraints, show them.
			else if( command.startsWith( "Show" ) ) {
				String alignments = backend.getDynalignAlignmentConstraints();
				handler.makeMessageDialog( alignments );
			}

			return;
		}

		// Determine the strand handled by this menu.
		String title = getText();
		str =
			( title.endsWith( "1" ) ) ? 1 :
			( title.endsWith( "2" ) ) ? 2 :
			0;

		// If the command comes from the "Chemical Modification" menu item,
		// choose a nucleotide to cleave.
		if( command.equals( "Chemical Modification" ) ) {
			new NucOrHelixDialog( command, "Base Number" ) {
				public String setConstraint( String data ) {
					Integer nuc = Integer.parseInt( data );
					return backend.setModifiedNucleotide( nuc, str );
				}
			};
		}

		// If the command comes from the maximum pairing distance menu item,
		// choose a new maximum pairing distance, if necessary.
		else if( command.contains( "Distance" ) ) {

			// Get the current max pairing distance and determine if it's set.
			int maxPair = backend.getMaxPair();
			String button = ( maxPair != -1 ) ? "Yes" : "No";
			if( maxPair == -1 ) { maxPair = 600; }

			// Create the input controls.
			IntegerField field = new IntegerField(
				"Maximum Distance", maxPair, 1 );
			String[] buttons = new String[]{ "Yes", "No" };
			RadioButtonPanel panel = RadioButtonPanel.makeHorizontal(
				"Limit Distance Between Paired Bases:", buttons );
			panel.setSelectedButton( button );

			// Put the input controls in a box and show the dialog.
			Box box = Box.createVerticalBox();
			box.add( panel );
			box.add( field );
			box.setPreferredSize(
				new Dimension( 300, box.getPreferredSize().height ) );
			String response = handler.makeTwoChoicesDialog( box );

			// If maximum pairing distance should be set, set that value.
			boolean setDistance =
				panel.getSelectedName().equals( "Yes" ) &&
				response.equals( "OK" );
			if( setDistance ) { backend.setMaxPair( field.getValue() ); }
		}

		// If the command comes from the "Double Stranded" menu item, choose a
		// nucleotide to force double stranded.
		else if( command.equals( "Double Stranded" ) ) {
			new NucOrHelixDialog( command, "Base Number" ) {
				public String setConstraint( String data ) {
					Integer nuc = Integer.parseInt( data );
					return backend.setDoubleStrandedNucleotide( nuc, str );
				}
			};
		}

		// If the command comes from the "FMN Cleavage" menu item, choose a
		// nucleotide to cleave.
		else if( command.equals( "FMN Cleavage" ) ) {
			new NucOrHelixDialog( command, "Base Number" ) {
				public String setConstraint( String data ) {
					Integer nuc = Integer.parseInt( data );
					return backend.setCleavedNucleotide( nuc, str );
				}
			};
		}

		// If the command comes from the "Force Pair" menu item, choose a
		// pair to force.
		else if( command.equals( "Force Base Pairs" ) ) {
			String helixHeaders = "Base 1;Base 2;Helix Length";
			new NucOrHelixDialog( command, helixHeaders ) {
				public String setConstraint( String data ) {
					String[] values = data.split( " " );
					Integer nuc1 = Integer.parseInt( values[0] );
					Integer nuc2 = Integer.parseInt( values[1] );
					Integer len = Integer.parseInt( values[2] );
					return backend.setForcedHelix( nuc1, nuc2, len, str );
				}
			};
		}

		// If the command comes from the "Maximum Loop" menu, choose a new
		// maximum loop size, if necessary.
		else if( command.contains( "Loop" ) ) {
			String fieldName =
				"<html>Maximum number of unpaired nucleotides<br/>" +
				"in internal/bulge loops:";
			int loop = backend.getMaxLoop();
			IntegerField field = new IntegerField( fieldName, loop, 0 );
			int newLoop = handler.makeInputDialog( field );
			backend.setMaxLoop( newLoop );
		}

		// If the command comes from the "Prohibit Base Pairs" menu item,
		// choose a pair to force.
		else if( command.equals( "Prohibit Base Pairs" ) ) {
			String helixHeaders = "Base 1;Base 2;Helix Length";
			new NucOrHelixDialog( command, helixHeaders ) {
				public String setConstraint( String data ) {
					String[] values = data.split( " " );
					Integer nuc1 = Integer.parseInt( values[0] );
					Integer nuc2 = Integer.parseInt( values[1] );
					Integer len = Integer.parseInt( values[2] );
					return backend.setProhibitedHelix( nuc1, nuc2, len, str );
				}
			};
		}

		// If the command is to reset constraints, do so.
		else if( command.startsWith( "Reset" ) ) {
			boolean clear = showDeletionConfirmation();
			if( clear ) { backend.clearFoldingConstraints( str ); }
		}

		// If the command is to read a constraints file, do so.
		else if( command.startsWith( "Restore" ) ) {
			String file = FileChooser.doOpen().getConstraints();
			String result = backend.readFoldingConstraintsFile( file, str );
			if( !result.equals( "" ) ) { handler.makeErrorDialog( result ); }
		}

		// If the command is to write a constraints file, do so.
		else if( command.startsWith( "Save" ) ) {
			String file = FileChooser.doSave().getConstraints();
			backend.writeFoldingConstraintsFile( file, str );
			handler.makeMessageDialog( "Folding constraints file written." );
		}

		// If the command comes from one of the SHAPE commands, show a SHAPE
		// constraints dialog.
		else if( command.contains( "SHAPE" ) ) {

			// Determine the SHAPE parameters.
			boolean isEnergy = command.contains( "Energy" );
			String desc1 = ( isEnergy ) ?
				"Slope (kcal/mol):" : "Threshold for Force Single Stranded:";
			double param1 = ( isEnergy ) ? 1.8 : 2;
			String desc2 = ( isEnergy ) ?
				"Intercept (kcal/mol):" :
				"Threshold for Chemical Modification:";
			double param2 = ( isEnergy ) ? -0.6 : 1;

			// Create the SHAPE file input field.
			final FileField file = FileField.createEnabled( "SHAPE Data File" );
			FilePanel filePanel = new FilePanel( new ActionListener() {
				public void actionPerformed( ActionEvent e ) {
					file.setText( FileChooser.doOpen().getSHAPE() );
				}
			}, file );
			filePanel.setPanelWidth( 400 );
			filePanel.makePanel();

			// Create the parameter values panel.
			DoubleField field1 = new DoubleField( desc1, param1 );
			DoubleField field2 = new DoubleField( desc2, param2 );
			NumberPanel values = new NumberPanel( field1, field2 );
			values.setPanelWidth( 400 );
			values.makePanel();

			// Put the input field and parameter values together into a panel.
			JPanel shapePanel = new JPanel( new BorderLayout() );
			shapePanel.add( filePanel, BorderLayout.NORTH );
			shapePanel.add( values, BorderLayout.SOUTH );

			// Show the SHAPE dialog. If constraints were selected, set them.
			String response = handler.makeTwoChoicesDialog( shapePanel );
			if( response.equals( "OK" ) ) {
				backend.setSHAPEFile( file.getText().trim() );
				backend.setSHAPEParam1( field1.getValue() );
				backend.setSHAPEParam2( field2.getValue() );
				backend.setSHAPEType( isEnergy );
			}
		}

		// If the command is to show constraints, show them.
		else if( command.startsWith( "Show" ) ) {
			String restraints = backend.getFoldingConstraints( str );
			handler.makeMessageDialog( restraints );
		}

		// If the command comes from the "Single Stranded" menu item, choose a
		// nucleotide to force single stranded.
		else if( command.equals( "Single Stranded" ) ) {
			new NucOrHelixDialog( command, "Base Number" ) {
				public String setConstraint( String data ) {
					Integer nuc = Integer.parseInt( data );
					return backend.setSingleStrandedNucleotide( nuc, str );
				}
			};
		}

		// If the command comes from the "Temperature" menu, choose a new
		// calculation temperature, if necessary.
		else if( command.contains( "Temperature" ) ) {
			double temperature = backend.getTemperature();
			DoubleField field =
				new DoubleField( "Temperature (degrees K):", temperature, 0 );
			double newTemperature = handler.makeInputDialog( field );
			backend.setTemperature( newTemperature );
		}
	}

	/**
	 * Show a dialog that asks whether constraints should be deleted.
	 *
	 * @return   True if constraints should be deleted, false if not.
	 */
	private boolean showDeletionConfirmation() {
		String msg = "This will reset all constraints.<br>Continue?";
		String reply = new SimpleDialogHandler().makeTwoChoicesDialog( msg );
		return reply.equals( "OK" );
	}

	/**
	 * An inner class which creates a dialog that sets constraints on single
	 * nucleotides or helices.
	 *
	 * @author Jessica S. Reuter
	 */
	public abstract class NucOrHelixDialog
		extends ValueSelectionDialog {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 *
		 * @param title    The title of the dialog.
		 * @param list     The list of field names, separated by semicolons.
		 */
		public NucOrHelixDialog( String title, String list ) {
			super( title );
			String[] titles = list.split( ";" );
			IntegerField[] numFields = new IntegerField[titles.length];
			for( int i = 1; i <= titles.length; i++ ) {
				numFields[i-1] = new IntegerField( titles[i-1], 1, 1,
					backend.getMaxConstraintIndex( str ) );
			}
			buildDialog( "Apply", numFields );
		}

		@Override
		public ActionListener createSelectionAction() {
			return new ActionListener() {
				public void actionPerformed( ActionEvent e ) {

					// Concatenate all the data fields into a string.
					String dataString = "";
					for( NumberField field: fields ) {
						dataString += ( field.getText() + " " );
					}

					// Set constraints, and show an error if one occurred.
					String result = setConstraint( dataString.trim() );
					if( !result.equals( "" ) ) {
						new SimpleDialogHandler().makeErrorDialog( result );
					}
				}
			};
		}

		/**
		 * Set a constraint using this dialog.
		 *
		 * @param data   The data string used to set the constraint.
		 */
		public abstract String setConstraint( String data );
	}
}
