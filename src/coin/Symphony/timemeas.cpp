/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2011 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#if !defined (_MSC_VER) && !defined (__DARWIN) && !defined (__MNO_CYGWIN) && !defined(__MINGW32__)
#include <sys/resource.h>
#endif
#include <stdio.h>

#include "sym_timemeas.h"

//extern int getrusage(int who, struct rusage *x);

double used_time(double *T)
{
#ifdef _OPENMP
   double t, oldT =  T ? *T : 0;
   struct timeval tp;

#if defined(_MSC_VER) || defined(__MNO_CYGWIN)
   clock_t tp_clock;
   long tp_long;
   tp_clock = clock();
   tp_long = (long) (tp_clock) / CLOCKS_PER_SEC;
   tp.tv_sec = tp_long;
   tp.tv_usec = 0;
#else
   (void) gettimeofday(&tp, NULL);
#endif

   t = (double)tp.tv_sec + ((double)tp.tv_usec)/1000000;
   if (T) *T = t;
   return (t - oldT);
#else

   /* FIXME: Windows CPU timing does not currently work */

#if defined(_MSC_VER) || defined (__MNO_CYGWIN) || defined(__MINGW32__)
   return (0);
#else

   double oldT =  *T;
   struct rusage x;

   (void) getrusage(RUSAGE_SELF, &x);
   *T = ((1e6 * (double) x.ru_utime.tv_sec) + (double)x.ru_utime.tv_usec);
   *T /= 1e6;
   return (*T - oldT);

#endif

/* END OF FIXME */

#endif
}

double wall_clock(double *T)
{
   double t, oldT =  T ? *T : 0;
   struct timeval tp;

#if defined(_MSC_VER) || defined (__MNO_CYGWIN)
   clock_t tp_clock;
   long tp_long;
   tp_clock = clock();
   tp_long = (long) (tp_clock) / CLOCKS_PER_SEC;
   tp.tv_sec = tp_long;
   tp.tv_usec = 0;
#else
   (void) gettimeofday(&tp, NULL);
#endif

   t = (double)tp.tv_sec + ((double)tp.tv_usec)/1000000;
   if (T) *T = t;
   return (t - oldT);
}

#if 0

#define MAX_SEC 100000000

void start_time(void)
{
   struct itimerval value = {{0, 0}, {MAX_SEC, 0}};
   struct itimerval ovalue = {{0, 0}, {0,0}};

   setitimer(ITIMER_VIRTUAL, &value, &ovalue);
}

double used_time(double *T)
{
   static struct itimerval value;

   getitimer(ITIMER_VIRTUAL, &value);
   return( ((double) MAX_SEC) -
	   ((double) value.it_value.tv_sec) -
	   ((double) value.it_value.tv_usec) / 1e6 - *T);
}

#endif
