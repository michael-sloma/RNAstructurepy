/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package RNAstructure_java_utilities;

import java.io.Serializable;

/**
 * A class that handles checking whether an OS is a Mac, and setting general
 * properties if that is the case.
 *
 * @author Jessica S. Reuter
 */
public class MacChecker
	implements Serializable {
	private static final long serialVersionUID = 20120802;

	/**
	 * Check the OS to see if it's a Mac, and if it is, set a few system
	 * properties to make things more "Mac-like."
	 */
	public static void checkOS() {
		if( notMac() ) { return; }
		System.setProperty( "apple.laf.useScreenMenuBar", "true" );
		System.setProperty( "com.apple.mrj.application.live-resize", "true" );
	}

	/**
	 * Get whether the current OS is a Mac.
	 *
	 * @return   True if the OS is a Mac, false if not.
	 */
	public static boolean isMac() {
		return ( System.getProperty( "mrj.version" ) != null );
	}

	/**
	 * Get whether the current OS is not a Mac.
	 *
	 * @return   True if the OS is not a Mac, false if it is.
	 */
	public static boolean notMac() {
		return ( isMac() == false );
	}
}
