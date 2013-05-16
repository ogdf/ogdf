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

#ifndef _BB_MACROS_H
#define _BB_MACROS_H

/*-------------------------- Random number generator ------------------------*/

#if defined(_MSC_VER) || defined (__MNO_CYGWIN) || defined(__MINGW32__) /* Different function call in
						  Windows */
#define SRANDOM(seed) srand(seed)
#define RANDOM() rand()
#else
#define SRANDOM(seed) srandom(seed)
#define RANDOM() random()
#endif

/*---------------------------- Allocation macros ----------------------------*/

#ifdef REMALLOC
#undef REMALLOC
#endif
#define REMALLOC(ptr, ptrtype, oldsize, newsize, block_size)	             \
{								             \
   if (!ptr || (oldsize < newsize)){				             \
      FREE(ptr);						             \
      oldsize = newsize + (int)(block_size);                                   \
      ptr = (ptrtype *) malloc((size_t)(oldsize) * sizeof(ptrtype));	     \
   }								             \
}

#ifdef REALLOC
#undef REALLOC
#endif
#define REALLOC(ptr, ptrtype, oldsize, newsize, block_size)		     \
{									     \
   if (!ptr || (oldsize < newsize)){					     \
      oldsize = newsize + (int)(block_size);				     \
      ptr = (ptrtype *) realloc((char *)ptr, (size_t)                        \
                                (oldsize * sizeof(ptrtype)));                \
   }									     \
}

/*---------------------------- PVM macros -----------------------------------*/

#define READ_INT_DESC(desc)						      \
{									      \
   receive_int_array(&(desc).size, 1);                                        \
   if ((desc).size > 0){						      \
      REMALLOC((desc).list, int, (desc).maxsize, (desc).size, BB_BUNCH);      \
      receive_int_array((desc).list, (desc).size);                            \
   }									      \
}

#define READ_CHAR_ARRAY_WITH_SIZE(cptr, cnum, maxcnum)	\
{								\
   receive_int_array(&cnum, 1);                                 \
   if (cnum > 0){						\
      REMALLOC(cptr, char, maxcnum, cnum, BB_BUNCH);		\
      receive_char_array(cptr, cnum);                           \
   }								\
}

#define READ_STR_LIST(snum, ssize, cptr, sptr)		\
{								\
   if (snum > 0){						\
      sptr = (char **) malloc(snum * sizeof(char *));		\
      cptr = (char *)  malloc(ssize * snum * CSIZE);		\
      receive_char_array(cptr, snum * ssize);              	\
      for (i = 0; i < snum; i++)				\
	 sptr[i] = cptr + i * ssize;				\
   }								\
}


/*--------------------- Parameter reading macros ----------------------------*/

#define READPAR_ERROR(x)						\
{									\
   (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", x);	\
   exit(1);								\
}

#define READ_INT_PAR(par)						\
if (sscanf(value, "%i", &(par)) != 1){					\
   (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", key);	\
   exit(1);								\
}

#define READ_STR_PAR(par)						\
if (sscanf(value, "%s", par) != 1){					\
   (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", key);	\
   exit(1);								\
}

#define READ_DBL_PAR(par)						\
if (sscanf(value, "%lf", &(par)) != 1){					\
   (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", key);	\
   exit(1);								\
}

#define READ_STRINT_PAR(par, str_array, array_size, value)		   \
{									   \
   for (i = array_size-1; i >= 0; i--){					   \
      if (! strcmp(str_array[i].str, value)){				   \
	 par |= str_array[i].code;					   \
	 break;								   \
      }									   \
   }									   \
   if (i < 0){								   \
      (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", key); \
      exit(1);								   \
   }									   \
}

/*------------------------ Copying macros -----------------------------------*/

#define COPY_DBL_ARRAY_DESC(newad, oldad)                                  \
if (newad.size > 0){                                                       \
   newad.stat = (int *) malloc(newad.size*ISIZE);                          \
   memcpy((char *)newad.stat, (char *)oldad.stat, oldad.size*ISIZE);       \
   if (newad.type == WRT_PARENT){                                          \
      newad.list = (int *) malloc(oldad.size*ISIZE);                       \
      memcpy((char *)newad.list, (char *)oldad.list, oldad.size*ISIZE);    \
   }                                                                       \
}

#define COPY_STAT(newad, oldad)                                            \
if (newad.size > 0){                                                       \
   newad.stat = (int *) malloc(newad.size*ISIZE);                          \
   memcpy((char *)newad.stat, (char *)oldad.stat, oldad.size*ISIZE);       \
}

#define COPY_ARRAY_DESC(newad, oldad)                                      \
newad = oldad;                                                             \
if (newad.size > 0){                                                       \
   newad.list = (int *) malloc(oldad.size*ISIZE);                          \
   memcpy((char *)newad.list, (char *)oldad.list, oldad.size*ISIZE);       \
}

/*--------------- Macro for calling user functions --------------------------*/

#define CALL_USER_FUNCTION(f)                                              \
switch (f){                                                                \
 case USER_ERROR:                                                          \
   printf("\n\n*********User error detected -- aborting***********\n\n");  \
   return(ERROR__USER);                                                    \
 default:                                                                  \
   break;                                                                  \
}

/*------------- Macro for calling wrapper functions -------------------------*/

#define CALL_WRAPPER_FUNCTION(f)                                           \
if ((termcode = f) < 0)                                                    \
   return(termcode);

/*---------------------- Standard macros ------------------------------------*/

#ifdef PRINT
#undef PRINT
#endif
#define PRINT(a, b, c) \
   if ((a) > (b)) printf c

#ifdef FREE
#undef FREE
#endif
#define FREE(p) if (p) {(void)free((char *)(p)); p = NULL;}

#ifndef MIN
#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))


#ifdef isset
#   undef isset
#   undef isclr
#endif
#ifdef ISSET_WITH_NBBY
#   ifndef NBBY
#      define NBBY 8
#   endif
#   define isset(a, i)     (((a)[(i)/NBBY] >> ((i) % NBBY)) & 1)
#else
#   define LOG_OF_BITS_PER_BYTE 3
#   define BITS_PER_BYTE_LESS_ONE 7
#   define isset(a, i)  (((a)[(i) >> LOG_OF_BITS_PER_BYTE] >> \
			  ((i) & BITS_PER_BYTE_LESS_ONE)) & 1)
#endif
#define isclr(a, i)  (! isset(a, i))

#ifdef setbit
#   undef setbit
#endif
#ifdef ISSET_WITH_NBBY
#   ifndef NBBY
#      define NBBY 8
#   endif
#   define setbit(a, i) ((a)[(i)/NBBY] |= (1 << ((i) % NBBY)))
#else
#   ifndef LOG_OF_BITS_PER_BYTE
#      define LOG_OF_BITS_PER_BYTE 3
#      define BITS_PER_BYTE_LESS_ONE 7
#   endif
#   define setbit(a, i) ((a)[(i) >> LOG_OF_BITS_PER_BYTE] |= \
			 (1 << ((i) & BITS_PER_BYTE_LESS_ONE)))
#endif

#endif
