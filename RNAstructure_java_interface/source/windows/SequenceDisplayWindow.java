/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package RNAstructure_java_interface.source.windows;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.BufferedInputStream;
import java.util.ArrayList;
import java.util.LinkedList;

import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.Clip;
import javax.sound.sampled.DataLine;
import javax.swing.Action;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.Timer;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.text.DefaultEditorKit;

import RNAstructure_java_interface.source.menus.FileMenu;
import RNAstructure_java_interface.source.menus.RolloverMenu;
import RNAstructure_java_interface.source.startup.ApplicationRootFrame;
import RNAstructure_java_utilities.BorderBuilder;
import RNAstructure_java_utilities.FileChooser;
import RNAstructure_java_utilities.ScrollerPane;
import RNAstructure_java_utilities.SimpleDialogHandler;

/**
 * A class responsible for displaying a window that allows a user to examine
 * raw sequence data.
 * <br><br>
 * Note that this class violates the usual 500 line class length coding standard
 * limit. This is because this window class contains menus unique to itself and
 * a variety of unique window actions. For the sake of organization and clarity
 * it is best to put all these together in one class.
 *
 * @author Jessica S. Reuter
 */
public class SequenceDisplayWindow
	extends InternalWindow {
	private static final long serialVersionUID = 20120802;

	/**
	 * The array of text areas in this window.
	 */
	private ArrayList<JTextArea> areas;

	/**
	 * A boolean, true if the sequence has been edited, false if not.
	 */
	private boolean edited = false;

	/**
	 * The position at which reading is occurring.
	 */
	private int position = 0;

	/**
	 * The button which can be clicked to start or stop reading.
	 */
	private JButton readButton;

	/**
	 * A boolean, true if reading while typing is active, false if not.
	 */
	private boolean reading = false;

	/**
	 * The timer that aids in reading nucleotides in rhythm.
	 */
	private Timer timer;

	/**
	 * Default Constructor.
	 * <br><br>
	 * This constructor is used when a blank display window is opened.
	 */
	public SequenceDisplayWindow() {

		// Set the title of this window and the main frame to reflect the
		// creation of a new sequence. Also, enable the menus and stop the
		// frame from closing automatically if the close button is clicked.
		setTitles( "New Sequence" );
		menuBar.enableMenus();
		setDefaultCloseOperation( SequenceDisplayWindow.DO_NOTHING_ON_CLOSE );

		// Insert the sequence display specific menu items in the file menu.
		FileMenu fileMenu = (FileMenu)menuBar.getMenu( 0 );
		fileMenu.insertMenuItem( "Save Sequence",
			"Save a sequence with its existing name.", this, 3 );
		fileMenu.insertMenuItem( "Save Sequence As...",
			"Save a sequence with a new name.", this, 4 );
		fileMenu.insertSeparator( 4 );

		// Initialize the text areas list.
		areas = new ArrayList<JTextArea>();

		// Create the input regions.
		BorderBuilder borders = new BorderBuilder();
		JPanel titlePane = buildInputRegion( 40, "Title:" );
		borders.makeEqualBorder( 5, titlePane );
		JPanel commentPane = buildInputRegion( 50, "Comment:" );
		borders.makeEqualBorder( 5, commentPane );
		JPanel sequencePane = buildInputRegion( 200,
			"Sequence: (Nucleotides in lower case are forced single " +
			"stranded in structure predictions.)" );
		borders.makeEqualBorder( 5, sequencePane );

		// Create the entire layout in a box and add the box to the window.
		Box box = Box.createVerticalBox();
		box.add( titlePane );
		box.add( commentPane );
		box.add( sequencePane );
		box.add( buildButtonPanel() );
		add( box );

		// Add a key listener to the sequence data input region so it restricts
		// nucleotides, and can read nucleotides out loud, if desired.
		getArea( 3 ).addKeyListener( new KeyAdapter() {

			/**
			 * Check to make sure a key is valid to be in the sequence region.
			 *
			 * @param e   The key event.
			 */
			public void checkKey( KeyEvent e ) {

				// Get the next character.
				char c = e.getKeyChar();

				// Check if the nucleotide is valid whitespace.
				boolean okSpace = Character.isWhitespace( c );

				// Check to see if the character is a valid nucleotide.
				boolean okNuc =
					c == 'a' || c == 'A' ||
					c == 'c' || c == 'C' ||
					c == 'g' || c == 'G' ||
					c == 't' || c == 'T' ||
					c == 'u' || c == 'U' ||
					c == 'x' || c == 'X';

				// Check to see if the character is a "movement" character.
				boolean okMove =
					( c == KeyEvent.VK_ENTER ) ||
					( c == KeyEvent.VK_BACK_SPACE );

				// If reading while typing is turned on, and the character is
				// a valid nucleotide, read it out.
				if( reading && okNuc ) { playClip( c ); }

				// If the character is not valid, don't output it.
				if( !( okSpace || okNuc || okMove ) ) { e.consume(); }
			}

			@Override
			public void keyPressed( KeyEvent e ) { checkKey( e ); }

			@Override
			public void keyTyped( KeyEvent e ) { checkKey( e ); }
		});

		// Add an internal frame listener that adds the following specialized
		// capabilities to this window:
		// 1. When the frame is activated, activate the save button.
		// 2. As the frame is closing, check if the sequence has been edited.
		// 3. When the frame is deactivated, deactivate the save button.
		final SequenceDisplayWindow saveActionHandler = this;
		final JButton save =
			(JButton)ApplicationRootFrame.getFrame().getToolBar()
			.getComponent( 3 );
		addInternalFrameListener( new InternalFrameAdapter() {

			@Override
			public void internalFrameActivated( InternalFrameEvent e ) {
				save.addActionListener( saveActionHandler );
				save.setEnabled( true );
			}

			@Override
			public void internalFrameClosing( InternalFrameEvent e ) {

				// If the window hasn't been edited, just dispose of the frame.
				if( !edited ) {
					e.getInternalFrame().dispose();
					return;
				}

				// Ask the user if they want to save the edited sequence.
				String msg =
					"The sequence has been modified.<br>Save Changes?";
				String response = dialogHandler.makeThreeChoicesDialog( msg );

				// If the sequence should be saved, save it.
				if( response.equals( "YES" ) ) {
					String file = saveSequence();
					if( !file.equals( "" ) ) { e.getInternalFrame().dispose(); }
				}

				// If the sequence shouldn't be saved, close the window.
				// Note that receiving a response of "CANCEL" does not close
				// the window; it only closes the dialog box opened above, as
				// no other action is attached to it.
				else if( response.equals( "NO" ) ) {
					e.getInternalFrame().dispose();
				}
			}

			@Override
			public void internalFrameDeactivated( InternalFrameEvent e ) {
				save.removeActionListener( saveActionHandler );
				save.setEnabled( false );
			}
		});

		// Initialize the reading timer.
		timer = new Timer( 750, this );
		timer.setActionCommand( "Timer" );
	}

	/**
	 * File Constructor.
	 * <br><br>
	 * This constructor is used when a sequence file is specified.
	 *
	 * @param file   The sequence file to open.
	 */
	public SequenceDisplayWindow( String file ) {

		// Call the blank window constructor to build the window before adding
		// the sequence.
		this();

		// Set the title of this window and the main frame to reflect the name
		// of the data file.
		setTitles( file );

		// Read in the data for the sequence, then set it in the window.
		backend.readSequenceData( file );
		getArea( 1 ).setText( backend.getSequenceTitle().trim() );
		getArea( 2 ).setText( backend.getSequenceComment().trim() );
		getArea( 3 ).setText( backend.getSequenceData().trim() );

		// Format the sequence and set it as unedited.
		formatSequence();
		edited = false;
	}

	/**
	 * Make the button panel at the bottom of the window.
	 * @return   The button panel.
	 */
	private JPanel buildButtonPanel() {

		// Create the button names array.
		String[] names = {
			"Format Sequence",
			"Read Sequence",
			"Fold as DNA",
			"Fold as RNA"
		};

		// Create the button panel, and save the second button as the reading
		// control button.
		JPanel buttonPanel = new JPanel( new GridLayout( 1, 0 ) );
		for( int i = 1; i <= 4; i++ ) {
			JButton button = new JButton( names[i-1] );
			button.addActionListener( this );
			buttonPanel.add( button );
			if( i == 2 ) { readButton = button; }
		}

		// Return the panel.
		return buttonPanel;
	}

	/**
	 * Build a text input region for this window. Each region is a panel
	 * consisting of a label and a text component.
	 *
	 * @param height   The height of the input region.
	 * @param text     The label text.
	 * @return         The input region panel.
	 */
	private JPanel buildInputRegion( int height, String text ) {

		// Create the input panel and add its description label.
		JPanel panel = new JPanel( new BorderLayout() );
		panel.add( new JLabel( text ), BorderLayout.NORTH );

		// Create the text area for this input region and set its defaults.
		JTextArea region = new JTextArea();
		region.setFont( new Font( "Monospaced", 0, 12 ) );
		region.setBackground( Color.WHITE );
		region.setWrapStyleWord( true );

		// Add a document listener to the text area so any change to the text
		// component sets the status as edited.
		region.getDocument().addDocumentListener( new DocumentListener() {
			public void changedUpdate( DocumentEvent e ) { edited = true; }
			public void insertUpdate( DocumentEvent e ) { edited = true; }
			public void removeUpdate( DocumentEvent e ) { edited = true; }
		});

		// Add the text area to a scroll pane and set its size, then add the
		// scroll pane to the panel and the text area to the array.
		ScrollerPane pane = new ScrollerPane( region, 100, height );
		panel.add( pane );
		areas.add( region );

		// Return the panel.
		return panel;
	}

	@Override
	protected void doActions( String command ) {

		// If the command comes from one of the "Fold" buttons, pop up a new
		// single strand folding window with the input already put in.
		if( command.startsWith( "Fold" ) ) {

			// Attempt to save the sequence, if edited, and if it wasn't saved,
			// return out of the possible folding.
			String file = getTitle();
			if( edited ) {
				file = saveSequence();
				if( file.equals( "" ) ) { return; }
			}

			// Check nucleic acid type, then show a window and set its data.
			dispose();
			String acid = ( command.endsWith( "RNA" ) ) ? "RNA" : "DNA";
			FoldSingleWindow window = new FoldSingleWindow( acid );
			window.viewWindow();
			window.setDataAutomatically( file );
		}

		// If the command comes from the "Format Sequence" button, format the
		// sequence into blocks so it's easier to read.
		else if( command.equals( "Format Sequence" ) ) { formatSequence(); }

		// If the command comes from the "Read Sequence" button, begin
		// reading the sequence.
		else if( command.equals( "Read Sequence" ) ) {
			position = 0;
			reading = true;
			timer.start();
			readButton.setText( "Stop Reading" );
		}

		// If the command is to save a sequence, save it either under its
		// existing name or under a new name.
		else if( command.startsWith( "Save" ) ) {
			boolean saveWithSameName =
				command.equals( "Save Sequence" ) &&
				!( getTitle().equals( "New Sequence" ) );
			if( saveWithSameName ) { saveSequence( getTitle() ); }
			else { saveSequence(); }
		}

		// If the command comes from the "Stop Reading" button, stop
		// reading the sequence.
		else if( command.equals( "Stop Reading" ) ) {
			reading = false;
			timer.stop();
			readButton.setText( "Read Sequence" );
		}

		// If the command comes from the timer, read a single nucleotide.
		else if( command.equals( "Timer" ) ) {
			try {

				// If the caret position is still within the length of the
				// sequence, highlight and read the next nucleotide.
				JTextArea area = getArea( 3 );
				if( position != area.getText().length() ) {

					// If the window was closed, stop the timer.
					if( !isVisible() ) { timer.stop(); }

					// Select the next nucleotide and play it if possible.
					area.getCaret().setSelectionVisible( false );
					area.select( position, position + 1 );
					Character base =
						area.getSelectedText().toUpperCase().charAt( 0 );
					if( base != ' ' ) {
						area.getCaret().setSelectionVisible( true );
						playClip( base );
					}
					position++;
				}

				// Otherwise if end has been reached, terminate the read.
				else {
					Thread.sleep( 600 );
					doActions( "Stop Reading" );
				}
			}

			// If a problem happened attempting to play a clip, show a
			// stack trace because the error won't be useful for the user.
			catch( Exception ex ) { ex.printStackTrace(); }
		}
	}

	/**
	 * Format the sequence in the text area.
	 */
	private void formatSequence() {

		// If reading while typing is on right now, turn it off.
		boolean isReading = reading;
		if( isReading ) { reading = false; }

		// Check if the sequence has already been edited.
		boolean wasEdited = edited;

		// Get unformatted text from the sequence area, as one long string with
		// all whitespace removed.
		JTextArea area = getArea( 3 );
		String text = area.getText().replaceAll( "[\\s]+", "" );
		int length = text.length();

		// Put blocks of the sequence 10 nucleotides long, separated by a
		// space, into the area, with a new line between every 5 blocks. If an
		// error occurs during this process, show it.
		try {

			// Create a list of lines 50 nucleotides long.
			LinkedList<String> lines = new LinkedList<String>();
			for( int i = 0; i < length; i += 50 ) {
				boolean fullLine = length - i >= 50;
				int end = ( fullLine ) ? i + 50 : i;
				String line = text.substring( i, end );
				if( fullLine ) { lines.add( line ); }
				else { lines.add( text.substring( i ) ); }
			}

			// For each line, split into blocks of 10 nucleotides, and append
			// each line onto one big formatted string.
			int size = lines.size();
			String fullText = "";
			for( int i = 0; i < size; i++ ) {
				StringBuilder builder = new StringBuilder( lines.get( i ) );
				for( int j = 10; j < builder.length(); j+=11 ) {
					builder.insert( j, ' ' );
				}
				fullText += ( builder.toString() + '\n' );
			}

			// Reset the pane text as the formatted text, and set the caret at
			// the beginning of the text.
			area.setText( fullText.trim() );
			area.setCaretPosition( 0 );
		} catch( Exception e ) {
			dialogHandler.makeErrorDialog( "Error formatting sequence." );
		}

		// Set the data of the window in the back end for writing later.
		backend.setSequenceTitle( getArea( 1 ).getText().trim() );
		backend.setSequenceComment( getArea( 2 ).getText().trim() );
		backend.setSequenceData( getArea( 3 ).getText().trim() );

		// If reading while typing was turned off before formatting the
		// sequence, turn it back on.
		if( isReading ) { reading = true; }

		// Set the edit state of the sequence to what it was before the
		// sequence was formatted. (That is, don't count formatting as a
		// sequence edit.)
		edited = wasEdited;
	}

	/**
	 * Get a text area from this window.
	 *
	 * @param index   The index of the area to retrieve, one-indexed.
	 * @return        The text area.
	 */
	private JTextArea getArea( int index ) { return areas.get( index - 1 ); }

	/**
	 * Play a sound clip, depending on the nucleotide.
	 *
	 * @param base   The nucleotide to play.
	 */
	private void playClip( char base ) {
		try {

			// Get the file name of the clip to play.
			String baseUp = Character.toString( base ).toUpperCase();
			String file =
				"/RNAstructure_java_interface/sounds/" + baseUp +  ".wav";

			// Get the audio information for the clip.
			BufferedInputStream buffer = new BufferedInputStream(
				SequenceDisplayWindow.class.getResourceAsStream( file ) );
			AudioInputStream input = AudioSystem.getAudioInputStream( buffer );
			DataLine.Info info =
				new DataLine.Info( Clip.class, input.getFormat() );

			// Play the clip.
			Clip clip = (Clip)AudioSystem.getLine( info );
			clip.open( input );
			clip.loop( 0 );
			Thread.sleep( 500 );
			clip.close();
		}

		// If a problem happened attempting to play a clip, show a stack
		// trace because the error won't be useful for the user.
		catch( Exception e ) {
			String message = "Error reading sequence out loud.";
			new SimpleDialogHandler().makeErrorDialog( message );
		}
	}

	/**
	 * Save a sequence.
	 *
	 * @return   The file that was saved.
	 */
	public String saveSequence() {
		String file = FileChooser.doSave().getSequence();
		if( file.equals( "" ) ) { return ""; }
		return saveSequence( file );
	}

	/**
	 * Save a sequence.
	 *
	 * @param file   The file name to save the sequence with.
	 * @return       The file that was saved.
	 */
	public String saveSequence( String file ) {

		// Format the sequence in the area.
		formatSequence();

		// Write the sequence file.
		if( file.endsWith( "fasta" ) ) {
			backend.writeFastaFile( file );
		} else {
			backend.writeSequenceFile( file );
		}

		// Set the titles of the window, set the state to be unedited, and
		// return the file.
		setTitles( file );
		edited = false;
		return file;
	}

	@Override
	protected RolloverMenu[] setMenus() {
		return new RolloverMenu[]{
			new EditMenu(), new ReadMenu()
		};
	}

	/**
	 * Set whether or not the window should be reading while typing.
	 *
	 * @param isReading   True if reading, false if not.
	 */
	public void setReading( boolean isReading ) { reading = isReading; }

	/**
	 * An inner class that creates a sequence editing menu.
	 *
	 * @author Jessica S. Reuter
	 */
	private class EditMenu
		extends RolloverMenu {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public EditMenu() {
			super( "Edit" );
			addMenuItem(
				"Cut", "Cut a block of text to the clipboard.", 'X' );
			addMenuItem(
				"Copy", "Copy a block of text to the clipboard.", 'C' );
			addMenuItem(
				"Paste", "Paste a block of text from the clipboard.", 'V' );
		}

		@Override
		protected void doMenuActions( String command ) {
			Action action =
				( command.equals( "Cut" ) ) ?
					new DefaultEditorKit.CutAction() :
				( command.equals( "Copy" ) ) ?
					new DefaultEditorKit.CopyAction() :
				new DefaultEditorKit.PasteAction();
			action.actionPerformed( null );
		}
	}

	/**
	 * An inner class that creates a sequence reading menu.
	 *
	 * @author Jessica S. Reuter
	 */
	private class ReadMenu
		extends RolloverMenu {
		private static final long serialVersionUID = 20120802;

		/**
		 * Constructor.
		 */
		public ReadMenu() {
			super( "Read" );
			addCheckBoxMenuItem( "Read While Typing",
				"Read a sequence out loud as it is typed into the " +
				"keyboard." );
		}

		@Override
		protected void doMenuActions( String command ) {
			boolean isReading = getItem( 0 ).isSelected();
			((SequenceDisplayWindow)ApplicationRootFrame.getFrame()
				.getMostRecentFrame()).setReading( isReading );
		}
	}
}
