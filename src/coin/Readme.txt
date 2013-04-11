***********************************************************
*             COIN-OR Library - Source Files              *
***********************************************************

The files in this directory are take from the COIN-OR project
(see http://www.coin-or.org) and are released under the
Eclipse Public License 1.0.

The sources are based on the following project releases:

- Clp 1.14.7
- Symphony 5.4.5

The files contain only minor modifications for splitting the
sources into header files (include/coin/) and implementation
files (src/coin) as well as a few compiler adaptions
(in particular for supporting MingGW on Windows).

Symphony/timemeas.cpp:
- conditional compilation for __MINGW32__ as for _MSC_VER

Symphony/master.cpp:
- replaced "(int)SYM_INFINITY" by INT_MAX to avoid compiler warning
- added include <limits.h>

Symphony/master_io.cpp:
- rename function "usage()" to "symusage()"
  (definition + a couple of usages, all only within this file)

Symphony/OsiSymSolverInterface.cpp:
- commented out some void code
- replaced some "(bool)value" casts by "value != 0"

CoinUtils/CoinModelUseful2.cpp:
- commented out some void code
