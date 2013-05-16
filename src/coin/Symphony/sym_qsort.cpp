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

#include "CoinSort.hpp"
#include <cstdlib>

/*===========================================================================*/
/*===========================================================================*/

static int cmpint(const void *p1, const void *p2){
   return (*((int *)p1) - *((int *)p2));
}

/*===========================================================================*/
/*===========================================================================*/

void qsort_i(int *bot, int nmemb)
{
  qsort (bot, nmemb, sizeof(int), cmpint);
}

/*===========================================================================*/
/*===========================================================================*/

void qsort_id(int *bot, double *botd, int nmemb)
{
   CoinSort_2(bot, bot+nmemb, botd);
}

/*===========================================================================*/
/*===========================================================================*/

void qsort_ic(int *bot, char *botc, int nmemb)
{
   CoinSort_2(bot, bot+nmemb, botc);
}

/*===========================================================================*/
/*===========================================================================*/

void qsort_ii(int *bot, int *bota, int nmemb)
{
   CoinSort_2(bot, bot+nmemb, bota);
}

/*===========================================================================*/
/*===========================================================================*/

void qsort_di(double *botd, int *boti, int nmemb)
{
   CoinSort_2(botd, botd+nmemb, boti);
}

/*===========================================================================*/
/*===========================================================================*/
/* calculate gcd of two integers i1, i2. */
/* TODO: replace with some function from CoinUtils */

int sym_gcd(int i1, int i2)
{
   int i;
   if (i1==0 && i2==0) {
      return 0;
   }

   if (i1<0) {
      i1 = -1*i1;
   }

   if (i2<0) {
      i2 = -1*i2;
   }

   if (i1==0) {
      return i2;
   }
   if (i2==0) {
      return i1;
   }

   while(1) {
      i = i2%i1;
      if (i==0) {
         return i1;
      } else {
         i2 = i1;
         i1 = i;
      }
   }
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/

