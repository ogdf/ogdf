***********************************************************
*             COIN-OR Library - Header Files              *
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

sym_macros.h:
- conditional compilation for __MINGW32__ as for _MSC_VER
