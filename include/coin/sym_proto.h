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

#ifndef _PROTO_H
#define _PROTO_H

#define MAX_FILE_NAME_LENGTH  80
#define MACH_NAME_LENGTH      80
#define MAX_LINE_LENGTH      255


#ifdef PROTO
#undef PROTO
#endif
/*#ifdef __GNUC__*/
#define PROTO(x) x
/*#else
#define PROTO(x) ()
#endif*/

#endif
