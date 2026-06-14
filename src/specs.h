#pragma once

/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2022 - 2024
*/

#ifdef _MSC_VER
#include ".\\gen\\toolbox.h"
#else
#include "./gen/toolbox.h"
#endif

// directory prefix
extern Str iox;		//        io directory prefix
extern Str bix;		// binary io directory prefix
extern Str tox;		//testoutput directory prefix
extern Str cdr;		// code root directory prefix

// directory
#ifdef _MSC_VER
const Str iopfx = "io\\";		//        io directory prefix
const Str bipfx = "bi\\";		// binary io directory prefix
const Str topfx = "tso\\";		//testoutput directory prefix
const Str dirslash = "\\";		// directory slash
const Str cdrpfx = "io\\..\\";	// code root directory prefix
#else
const Str iopfx = "../io/";		//        io directory prefix
const Str bipfx = "bi/";		// binary io directory prefix
const Str topfx = "tso/";		//testoutput directory prefix
const Str cdrpfx = "io/../";	// code root directory prefix
const Str dirslash = "/";		// directory slash
#endif
inline void iox_exist()
{
	if (!dir_exist(iox)) ERR(NAV(iox) + " does not exist!");
}
inline void io_init()
{
	if (!dir_exist(iox)) {
		WRN(NAV(iox) + " does not exist!");
		int flag = MKDIR(iox);
		if (flag == 0) {
			WRN(NAV(iox) + " is created.");
		}
		else {
			WRN(NAV(iox) + "creation fails.");
		}
	}
	
	if (!dir_exist(bix)) {
		WRN(NAV(bix) + " does not exist!");
		int flag = MKDIR(bix);
		if (flag == 0) {
			WRN(NAV(bix) + " is created.");
		}
		else {
			WRN(NAV(bix) + "creation fails.");
		}
	}

	if (!dir_exist(tox)) {
		WRN(NAV(tox) + " does not exist!");
		int flag = MKDIR(tox);
		if (flag == 0) {
			WRN(NAV(tox) + " is created.");
		}
		else {
			WRN(NAV(tox) + "creation fails.");
		}
	}
}

