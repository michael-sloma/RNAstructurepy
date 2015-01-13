/*
 * A header file for a progress dialog used with the RNAstructure GUI.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef SWIG_TPROGRESSDIALOG_H
#define SWIG_TPROGRESSDIALOG_H

#include "ProgressMonitor.h"

class TProgressDialog {
 public:
  // Public constructor and method

  /*
   * Name:        Constructor
   * Description: Constructs the TProgressDialog proxy object.
   * Argument:
   *          1. A reference to a ProgressMonitor
   */
  TProgressDialog( ProgressMonitor& monitor );

  /*
   * Name:        update
   * Description: Updates the underlying C++ class.
   * Argument:
   *     1.   The new percent value the underlying C++ class should hold.
   */
  void update( int percent );

 private:
  // The progress monitor
  ProgressMonitor* progressMonitor;
};

#endif /* SWIG_TPROGRESS_DIALOG_H */
