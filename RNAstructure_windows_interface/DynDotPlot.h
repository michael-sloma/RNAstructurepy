#pragma once
#include "plotdoc.h"

class CDynDotPlot :
	public PlotDoc
{
	DECLARE_DYNCREATE(CDynDotPlot)
public:
	CDynDotPlot(void);
	CDynDotPlot(CString Filename, int seq, CRNAstructureApp *app);
	~CDynDotPlot(void);
};
