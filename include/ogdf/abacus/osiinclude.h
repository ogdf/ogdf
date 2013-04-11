#ifndef OSI_INCLUDE_H
#define OSI_INCLUDE_H

#ifdef OSI_CBC
#include <coin/OsiCbcSolverInterface.hpp>
#endif
#ifdef OSI_CLP
#include <coin/OsiClpSolverInterface.hpp>
#endif
#ifdef OSI_CPX
#include <coin/OsiCpxSolverInterface.hpp>
#endif
#ifdef OSI_DYLP
#include <coin/OsiDylpSolverInterface.hpp>
#endif
#ifdef OSI_FORTMP
#include <coin/OsiFmpSolverInterface.hpp>
#endif
#ifdef OSI_GLPK
#include <coin/OsiGlpkSolverInterface.hpp>
#endif
#ifdef OSI_MOSEK
#include <coin/OsiMskSolverInterface.hpp>
#endif
#ifdef OSI_OSL
#include <coin/OsiOslSolverInterface.hpp>
#endif
#ifdef OSI_SOPLEX
#include <coin/OsiSpxSolverInterface.hpp>
#endif
#ifdef OSI_SYM
#include <coin/OsiSymSolverInterface.hpp>
#endif
#ifdef OSI_VOL
#include <coin/OsiVolSolverInterface.hpp>
#endif
#ifdef OSI_XPRESS
#include <coin/OsiXprSolverInterface.hpp>
#endif
#ifdef OSI_GRB
#include <coin/OsiGrbSolverInterface.hpp>
#endif
#ifdef OSI_CSDP
#include <coin/OsiCsdpSolverInterface.hpp>
#endif

#include <coin/OsiSolverInterface.hpp>

#endif // end ifndef
