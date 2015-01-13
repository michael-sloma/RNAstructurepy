/*
 * (c) 2011 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JLabel;

/**
 * A class that creates images in various contexts for the RNAstructure GUI.
 *
 * @author Jessica S. Reuter
 */
public class ImageGrabber
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Create an image.
	 *
	 * @param path           The path to the image.
	 * @return               The image.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static BufferedImage getImage( String path )
		throws IOException {
		InputStream iconIn = ImageGrabber.class.getResourceAsStream( path );
		return ImageIO.read( iconIn );
	}

	/**
	 * Create an image icon.
	 *
	 * @param path           The path to the image.
	 * @return               The image icon.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static ImageIcon getImageIcon( String path )
		throws IOException {
		return new ImageIcon( getImage( path ) );
	}

	/**
	 * Create an image label.
	 *
	 * @param path           The path to the image.
	 * @return               The image label.
	 * @throws IOException   If an image cannot be read correctly. Since all
	 *                       image paths are hard coded, this should never
	 *                       be thrown.
	 */
	public static JLabel getImageLabel( String path )
		throws IOException {
		return new JLabel( getImageIcon( path ) );
	}
}
