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

#ifndef _QSORT_H
#define _QSORT_H

//void qsort(char *bot, unsigned int nmemb, int size,
//	   int (*compar)(const void *, const void *));

void qsort_i(int *bot, int nmemb);
void qsort_id(int *bot, double *botd, int nmemb);
void qsort_ic(int *bot, char *botc, int nmemb);
void qsort_ii(int *bot, int *bota, int nmemb);
void qsort_di(double *botd, int *boti, int nmemb);
/* TODO: replace with some function from CoinUtils */
int sym_gcd(int i1, int i2);

#endif
