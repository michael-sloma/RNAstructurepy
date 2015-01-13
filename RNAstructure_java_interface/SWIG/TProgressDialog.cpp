/*
 * An implementation file for a progress dialog used with the RNAstructure GUI.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "TProgressDialog.h"

// Constructor
TProgressDialog::TProgressDialog( ProgressMonitor& monitor ) {
  progressMonitor = &monitor;
}

// Update the value of the monitor.
void TProgressDialog::update( int percent ) {
  progressMonitor->setMonitorProgress( percent );
}
