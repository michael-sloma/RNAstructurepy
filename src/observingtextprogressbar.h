/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef OBSERVINGTEXTPROGRESSBAR_H
#define OBSERVINGTEXTPROGRESSBAR_H

#ifdef _JAVA_GUI
#include "../RNAstructure_java_interface/SWIG/TProgressDialog.h"
#else

#ifndef _WINDOWS
#include "TProgressDialog.h"
#endif // !_WINDOWS

#endif // JAVA GUI                                                                                                                                                                                        

#include "observer.h"

class ObservingTextProgressBar : public TProgressDialog,
                                 public Observer {
private:
  int value;
  int max;
  
public:
  ObservingTextProgressBar(int _max = 100, bool Silent = false);
  void setMaximum(int _max);
  void notify();
  bool silent;
};

#endif
