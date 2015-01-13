/*
 * A header file for a progress monitor used with the RNAstructure GUI.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef PROGRESSMONITOR_H
#define PROGRESSMONITOR_H

class ProgressMonitor {
 public:
  // Public constructor and methods

  /*
   * Name:        Constructor
   * Description: Initializes the integer holding progress.
   */
  ProgressMonitor();

  /*
   * Name:        getMonitorProgress
   * Description: Get how far along an underlying task has progressed.
   * Return:
   *     An integer which is the percent progress, from 0 to 100.
   */
  int getMonitorProgress();

  /*
   * Name:        setMonitorProgress
   * Description: Set how far along an underlying task has progressed.
   * Argument:
   *          1. The new percent progress value to set
   */
  void setMonitorProgress( int percent );

 private:
  // Private variable: percent of total progress completed
  int percentProgress;
};

#endif /* PROGRESSMONITOR_H */
