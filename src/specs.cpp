/*
coded by Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#include "specs.h"

// directory
#ifdef _MSC_VER
extern Str iox = iopfx;		//        io directory prefix
extern Str bix = bipfx;		// binary io directory prefix
extern Str tox = topfx;		// binary io directory prefix
extern Str cdr = cdrpfx;	// code root directory prefix
#else
extern Str iox = iopfx;		//        io directory prefix
extern Str bix = bipfx;		// binary io directory prefix
extern Str tox = topfx;		// binary io directory prefix
extern Str cdr = cdrpfx;	// code root directory prefix
#endif
