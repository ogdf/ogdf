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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sym_constants.h"
#include "sym_primal_heuristics.h"
#include "sym_macros.h"

/* Functions related to solution pool */
#ifdef COMPILE_IN_LP

/*===========================================================================*/
/*===========================================================================*/
int sp_add_solution (lp_prob *p, int cnt, int *indices, double *values,
      double obj_value, int bc_index)
{
   sp_desc *sp = p->tm->sp;
   sp_solution *sol;

   //TODO: check duplicates

   if (sp->num_solutions == sp->max_solutions &&
         sp->solutions[0]->objval>=obj_value+p->lp_data->lpetol) {
      /* delete first solution and move everything up by 1 */
      sp_delete_solution(sp,0);
      /*
      for (i=0;i<(sp->max_solutions-1);i++) {
         sp->solutions[i] = sp->solutions[i+1];
      }
      */
   } else if (sp->num_solutions == sp->max_solutions) {
      /* pool is full and the new solution is worse than any stored solution
       */
      return 0;
   }
   sol = sp->solutions[sp->num_solutions];
   sol->objval = obj_value;
   sol->xlength = cnt;
   sol->xind = (int *) malloc(ISIZE*cnt);
   memcpy(sol->xind,indices,ISIZE*cnt);
   sol->xval = (double *) malloc(DSIZE*cnt);
   memcpy(sol->xval,values,DSIZE*cnt);
   sol->node_index = bc_index;
   sp->num_solutions++;
   sp->total_num_sols_found++;
   PRINT(p->par.verbosity,5,("sp: solution pool size = %d \n",
            sp->num_solutions));
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_delete_solution (sp_desc *sp, int position)
{
   int i;
   if (position>=sp->num_solutions) {
      return 0;
   }

   FREE(sp->solutions[position]->xind);
   FREE(sp->solutions[position]->xval);
   for (i=position; i<sp->num_solutions-1; i++) {
      sp->solutions[i]->xind=sp->solutions[i+1]->xind;
      sp->solutions[i]->xval=sp->solutions[i+1]->xval;
      sp->solutions[i]->objval = sp->solutions[i+1]->objval;
      sp->solutions[i]->xlength = sp->solutions[i+1]->xlength;
      sp->solutions[i]->node_index = sp->solutions[i+1]->node_index;
   }
   sp->solutions[sp->num_solutions-1]->xlength = 0;
   sp->num_solutions--;
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_is_solution_in_sp (lp_prob *p, int cnt, int *indices, double *values,
      double obj_value)
{
   /* not implemented yet */
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_initialize(tm_prob *tm)
{
   int i;
   tm->sp = (sp_desc*)malloc(sizeof(sp_desc));
   sp_desc *sp = tm->sp;
   sp->max_solutions = 10;
   sp->num_solutions = 0;
   sp->total_num_sols_found = 0;
   sp->solutions = (sp_solution **) malloc(sp->max_solutions*sizeof(sp_solution*));
   for (i=0;i<sp->max_solutions;i++) {
      sp->solutions[i] = (sp_solution *) malloc(sizeof(sp_solution));
   }

   /* TODO: put the above as parameters */

   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_free_sp(sp_desc *sp)
{
   int i;

   //TODO: move this to master_io()
   //printf("Total number of solutions found = %d\n",sp->total_num_sols_found);
   for (i=sp->num_solutions-1; i>=0; i--) {
      sp_delete_solution(sp,i);
   }
   for (i=sp->max_solutions-1; i>-1; i--) {
      FREE(sp->solutions[i]);
   }
   FREE(sp->solutions);
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
#endif //COMPILE_IN_LP
