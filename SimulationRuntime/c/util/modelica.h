/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */

/*! \file modelica.h
 * Description: This is the C header file for the C code generated from
 * Modelica. It includes e.g. the C object representation of the builtin types
 * and arrays, etc.
 */

#ifndef MODELICA_H_
#define MODELICA_H_

#if defined(__cplusplus)
extern "C" {
#endif

/* adrpo: extreme windows crap! */
#if defined(__MINGW32__) || defined(_MSC_VER)
#define DLLImport   __declspec( dllimport )
#define DLLExport   __declspec( dllexport )
#else
#define DLLImport /* extern */
#define DLLExport /* nothing */
#endif

#if defined(IMPORT_INTO)
#define DLLDirection DLLImport
#else /* we export from the dll */
#define DLLDirection DLLExport
#endif

#if defined(__MINGW32__) || defined(_MSC_VER)
 #define WIN32_LEAN_AND_MEAN
#if !defined(NOMINMAX)
 #define NOMINMAX
 #include <windows.h>
#endif
#endif

#include <stdlib.h>
#include <limits.h>

#if defined(__MINGW32__) || defined(_MSC_VER)
#define EXIT(code) exit(code)
#else
/* We need to patch exit() on Unix systems
 * It does not change the exit code of simulations for some reason! */
#include <unistd.h>
#define EXIT(code) {fflush(NULL); _exit(code);}
#endif

#include "omc_inline.h"

#include "util/modelica_string.h"
#include "util/memory_pool.h"
#include "util/index_spec.h"

#include "util/string_array.h"
#include "util/boolean_array.h"

#include "util/real_array.h"
#include "util/integer_array.h"
#include "util/generic_array.h"

#include "util/utility.h"
#include "util/division.h"

#include <assert.h>
#include "util/read_write.h"
#include "simulation_data.h"
#include "meta/meta_modelica.h"
#include "meta/meta_modelica_builtin.h"
#include "util/varinfo.h"


/* math functions (-lm)*/

/* Special Modelica builtin functions*/
#define smooth(P,EXP)    (EXP)
#define semiLinear(x,positiveSlope,negativeSlope) (x>=0?positiveSlope*x:negativeSlope*x)

/* sign function */
#define sign(v) (v>0?1:(v<0?-1:0))

#if defined(_MSC_VER)
#define fmax(x, y) ((x>y)?x:y)
#define fmin(x, y) ((x<y)?x:y)
#define snprintf sprintf_s
#endif

#if defined(__cplusplus)
}
#endif

#endif
