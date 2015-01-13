/*
 * An implementation file for a progress monitor used with the RNAstructure GUI.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#include "ProgressMonitor.h"

// Constructor
ProgressMonitor::ProgressMonitor() { percentProgress = 0; }

// Accessor
int ProgressMonitor::getMonitorProgress() { return percentProgress; }

// Mutator
void ProgressMonitor::setMonitorProgress( int percent ) {
  percentProgress = percent;
}
