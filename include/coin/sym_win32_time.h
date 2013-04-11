#ifndef WIN32_TIME_H
#define WIN32_TIME_H

#if defined(_MSC_VER) || defined (__MNO_CYGWIN)
#include <time.h>

struct timeval {
    long    tv_sec;         /* seconds */
    long    tv_usec;        /* and microseconds */
};

struct rusage {
       struct timeval   ru_utime;   /* user time used */
       struct timeval   ru_stime;   /* system time used */
       long             ru_maxrss;  /* maximum resident set size */
       long             ru_idrss;   /* integral resident set size */
       long             ru_minflt;  /* page faults not requiring physical I/O*/
       long             ru_majflt;  /* page faults requiring physical I/O */
       long             ru_nswap;   /* swaps */
       long             ru_inblock; /* block input operations */
       long             ru_oublock; /* block output operations */
       long             ru_msgsnd;  /* messages sent */
       long             ru_msgrcv;  /* messages received */
       long             ru_nsignals;/* signals received */
       long             ru_nvcsw;   /* voluntary context switches */
       long             ru_nivcsw;  /* involuntary context switches */
};

#pragma warning(disable:4244)
#pragma warning(disable:4013)
#pragma warning(disable:4068)
#pragma warning(disable:4761)
#pragma warning(disable:4305)
#pragma warning(disable:4101)
#pragma warning(disable:4700)
#pragma warning(disable:4018)


#endif
#endif
