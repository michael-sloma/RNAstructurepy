/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "observingtextprogressbar.h"

#include "observer.h"

#ifdef _JAVA_GUI
#include "../RNAstructure_java_interface/SWIG/TProgressDialog.h"
#else

#ifndef _WINDOWS
#include "TProgressDialog.h"
#endif // !_WINDOWS

#endif // JAVA GUI

ObservingTextProgressBar::ObservingTextProgressBar(int _max, bool Silent)
  : value(0),
    max(_max),
	silent(Silent)
{
}

void ObservingTextProgressBar::setMaximum(int _max) {
  max = _max;
}

void ObservingTextProgressBar::notify() {
	if (!silent) update((100 * ++value) / max);
}
