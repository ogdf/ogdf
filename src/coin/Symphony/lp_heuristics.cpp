/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Ashutosh Mahajan                               */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sym_qsort.h"
#include "sym_lp.h"
#include "sym_constants.h"
#include "sym_lp_solver.h"
#include "sym_primal_heuristics.h"
#include "sym_macros.h"

/*===========================================================================*/
/*===========================================================================*\
 * This file contains heuristics to find an integral solution after an LP
 * is solved.
 \*===========================================================================*/
/*===========================================================================*/
/*
 * TODO:
 * make independent of solver
 */

int feasibility_pump (lp_prob *p, char *found_better_solution,
      double &solution_value, double *betterSolution)
{
   int                      termcode    = FUNCTION_TERMINATED_NORMALLY;
   FPdata                  *fp_data     = (FPdata*) malloc(sizeof(FPdata));
   LPdata                  *lp_data     = p->lp_data;
   LPdata                  *new_lp_data = (LPdata *)calloc(1,sizeof(LPdata));
   /* no. of max pumping cycles */
   int                      max_iter    = p->par.fp_max_cycles;
   int                      n           = lp_data->n;
   /* use OSI to get lp data */
   OsiSolverInterface      *model       = p->lp_data->si;
   const CoinPackedMatrix  *matrix      = model->getMatrixByRow();
   const double            *lp_r_low    = model->getRowLower();
   const double            *lp_r_up     = model->getRowUpper();
   int                      i, r, iter, cnt, verbosity;
   int                     *indices;
   double                  *values;
   double                   fp_time, last_fp_time, real_obj_value, target_ub;
   FPvars                 **vars;
   double                   gap           = model->getInfinity();
   double                   obj_lb        = lp_data->objval;
   double                   total_time    = 0;
   const double            *mip_obj       = model->getObjCoefficients();
   char                     is_feasible   = FALSE;
   double                  *x_ip, *x_lp, new_solution_value;
   const double             fp_display_interval = p->par.fp_display_interval;
   /* number of solutions with obj value more than the best */
   int                      num_poor_sols = 0;
   int                      num_better_sols = 0;
   const double             lpetol = p->lp_data->lpetol;
   int                      fp_poor_sol_lim = p->par.fp_poor_sol_lim_fac;
   int                      total_iter_cnt = 0;
   fp_time                                = used_time(&total_time);

   if (p->lp_stat.fp_calls < 1) {
      CoinSeedRandom(17000);
   }

   /* total_time and fp_time both now have total time used by symphony's lp
    * process */
   fp_time                                = used_time(&total_time);
   last_fp_time                           = fp_time;
   /* fp_time should now be zero and total_time be still the same */

   verbosity = fp_data->verbosity         = p->par.verbosity;
   if (p->bc_index<1) {
      PRINT(verbosity, 0, ("starting feasibility pump\n"));
   }

   *found_better_solution = FALSE;
   fp_data->mip_obj       = (double *)malloc(n*DSIZE);
   fp_data->flip_fraction = p->par.fp_flip_fraction;
   fp_data->sos_row_filled = 0;
   fp_data->sos_var_fixed_zero = 0;
   fp_data->can_check_sos = FALSE;

   if(p->mip->matbeg && p->mip->mip_inf &&
      p->mip->mip_inf->binary_sos_row_num > 0){
      fp_data->can_check_sos = TRUE;
      fp_data->sos_row_filled = (char *)malloc(p->mip->m*CSIZE);
      //fp_data->sos_var_fixed_zero = (char *)malloc(p->mip->n*CSIZE);
   }

   memcpy(fp_data->mip_obj,mip_obj,n*DSIZE);

   /* initialize the lp solver. load the current basis */
   fp_initialize_lp_solver(p, new_lp_data, fp_data);
   x_ip = fp_data->x_ip;
   x_lp = fp_data->x_lp;

   if (p->has_ub) {
      solution_value = p->ub-p->par.granularity;
   }
   else {
      solution_value = model->getInfinity();
   }

   if (p->has_ub && p->mip->mip_inf &&
         (p->mip->mip_inf->obj_size <= p->mip->mip_inf->max_row_size ||
          p->mip->mip_inf->obj_size < n/10)) {
      fp_add_obj_row(new_lp_data, n, mip_obj, p->ub-p->par.granularity);
   }
   /* round the x_lp and store as x_ip, it will usually become infeasible */
   vars = fp_data->fp_vars;

   /* do the following max_iter times */
   fp_time += used_time(&total_time);
   int fp_override_cnt = 0;
   /*
   if(p->lp_stat.fp_calls == 1){
      p->par.fp_time_limit += 10;
      p->par.fp_max_initial_time += 10;
   }else if (p->lp_stat.fp_calls == 2){
      p->par.fp_time_limit -= 10;
      p->par.fp_max_initial_time -= 10;
   }
   */
   if(p->lp_stat.fp_calls < 1){
      p->par.fp_time_limit += 20;
   }else if(p->lp_stat.fp_calls < 2){
      p->par.fp_time_limit -= 20;
      p->par.fp_max_initial_time += 20;
   }else if(p->lp_stat.fp_calls < 3){
      p->par.fp_max_initial_time -= 20;
   }

   for (iter=0; (iter<max_iter && fp_time<p->par.fp_time_limit &&
		 fp_time + p->comp_times.fp < p->par.fp_max_initial_time) ||
	   fp_override_cnt > 0; iter++) {
      if (fp_time - last_fp_time > fp_display_interval || verbosity > 5) {
         PRINT(verbosity, 0,
               ("feasibility pump: starting iteration %d, time used = %.2f\n",
                iter, fp_time));
         last_fp_time = fp_time;
      }

      is_feasible = FALSE;
      /* solve an lp */
       fp_round(p, fp_data, new_lp_data);
      if (fp_data->x_bar_len[fp_data->iter] == -1) {
         /*
          * the cost and reference point are same as some other iteration. we
          * should stop here because we are cycling
          */
         PRINT(verbosity,5,("fp: leaving because of cycling\n"));
         fp_data->iter++;
         break;
      }
      fp_is_feasible (lp_data,matrix,lp_r_low,lp_r_up,fp_data,&is_feasible);

      if (is_feasible == TRUE) {
         new_solution_value = 0;
         for (i=0;i<n;i++) {
            new_solution_value += x_ip[i]*mip_obj[i];
         }
         if (new_solution_value<solution_value-p->par.granularity-lpetol) {
	    /* we found what we wanted */
	    memcpy(betterSolution, x_ip, n*DSIZE);

	    /* we found what we wanted */
	    memcpy(betterSolution, x_ip, n*DSIZE);

	    solution_value = new_solution_value;
            indices = p->lp_data->tmp.i1;          /* n */
            values  = p->lp_data->tmp.d;           /* n */
            cnt     = collect_nonzeros(p, betterSolution, indices, values);
            gap     = (solution_value -
                      obj_lb)/(fabs(solution_value)+0.001)*100;
            p->lp_stat.fp_num_sols++;
            num_better_sols++;
            PRINT(verbosity,5,("fp: found solution with value = %f\n",
                     solution_value));
            PRINT(verbosity,5,("fp: gap = %f\n", gap));
            sp_add_solution(p,cnt,indices,values,
                  solution_value+p->mip->obj_offset,p->bc_index);
            if (gap <= p->par.fp_min_gap) {
               *found_better_solution = TRUE;
               fp_data->iter++;
               break;
            }
            target_ub = (obj_lb + solution_value)/2;
            if (p->mip->mip_inf && (p->mip->mip_inf->obj_size <=
                     p->mip->mip_inf->max_row_size
                  || p->mip->mip_inf->obj_size < n/10)) {
               if (*found_better_solution != TRUE && p->has_ub==FALSE) {
                  // add another objective function constraint to lower the
                  // objective value.
                  fp_add_obj_row(new_lp_data, n, mip_obj, target_ub);
               } else {
                  r = new_lp_data->m-1;
                  change_rhs(new_lp_data, 1, &r, &target_ub);
               }
            }
            *found_better_solution = TRUE;
            fp_poor_sol_lim = p->par.fp_poor_sol_lim_fac *
                              num_better_sols;
	    /* menal ---*/
	    if(p->bc_level > 0) {
              fp_data->iter++;
              break;
            }
	    /* --- */
         } else {
            num_poor_sols++;
            /*
            PRINT(verbosity,5,("fp: rejecting poor solution with value = %f\n",
                     solution_value));
            PRINT(verbosity,5,("fp: number of poor sols = %d, better sols = %d, limit=%d\n",
                     num_poor_sols, num_better_sols, fp_poor_sol_lim));
            */
            if (num_poor_sols > fp_poor_sol_lim) {
            /*
               PRINT(verbosity,5,("fp: breaking because of too many (%d) poor"
                       " solutions\n", num_poor_sols));
            */
               fp_data->iter++;
               break;
            }
         }
      }

      PRINT(verbosity,5,("fp: solve lp %d\n",iter));
      p->lp_stat.lp_calls++;
      p->lp_stat.fp_lp_calls++;

      if (fp_solve_lp(new_lp_data, fp_data, &is_feasible) !=
            FUNCTION_TERMINATED_NORMALLY) {
         break;
      }

      fp_data->iter++;
      fp_time += used_time(&total_time);
      total_iter_cnt += fp_data->iterd;

   }

   p->lp_stat.fp_poor_sols = num_poor_sols;
   p->lp_stat.fp_lp_total_iter_num += total_iter_cnt;
   close_lp_solver(new_lp_data);
   /* free all the allocated memory */
   FREE(new_lp_data->x);
   FREE(new_lp_data->lb);
   FREE(new_lp_data->ub);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->dualsol);
   FREE(new_lp_data->dj);
   FREE(new_lp_data->tmp.c);
   FREE(new_lp_data->tmp.d);
   FREE(new_lp_data->tmp.i1);
   FREE(new_lp_data);

   for (i=0;i<n;i++) {
      FREE(fp_data->fp_vars[i]);
   }

   for (i=0;i<fp_data->iter;i++) {
      FREE(fp_data->x_bar_val[i]);
      FREE(fp_data->x_bar_ind[i]);
   }
   FREE(fp_data->x_bar_val);
   FREE(fp_data->x_bar_ind);
   FREE(fp_data->x_bar_len);
   FREE(fp_data->fp_vars);
   FREE(fp_data->sos_row_filled);
   FREE(fp_data->sos_var_fixed_zero);
   FREE(fp_data->obj);
   FREE(fp_data->mip_obj);
   FREE(fp_data->x_lp);
   FREE(fp_data->x_ip);
   FREE(fp_data->index_list);
   FREE(fp_data->alpha_p);
   FREE(fp_data);

   /* update stats */
   fp_time                        += used_time(&total_time);
   p->comp_times.fp               += fp_time;
   p->lp_stat.fp_calls++;
   if (*found_better_solution==TRUE) {
      if (p->mip->obj_sense == SYM_MAXIMIZE){
         real_obj_value=-solution_value+p->mip->obj_offset;
      } else {
         real_obj_value=solution_value+p->mip->obj_offset;
      }
      PRINT(verbosity,5,("fp: found solution = %10.2f time = %10.2f\n",
               real_obj_value,total_time));
   }

   if (p->bc_index<1 || verbosity > 5) {
      PRINT(verbosity, 0, ("leaving feasibility pump.\n"));
   }

   return termcode;
}


/*===========================================================================*/
int fp_is_feasible (LPdata *lp_data, const CoinPackedMatrix *matrix,
		    const double *r_low, const double *r_up, FPdata *fp_data,
		    char *is_feasible )
{
   /* check if x is a integer feasible solution to problem in p */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double lpetol = lp_data->lpetol;
   //int n = fp_data->n0;
   int m = fp_data->m0;
   //FPvars **vars = fp_data->fp_vars;
   int i,c,j;
   double Ractivity;
   const int *r_matbeg = matrix->getVectorStarts();
   const int *r_matlen = matrix->getVectorLengths();
   const int *r_matind = matrix->getIndices();
   const double *r_matval = matrix->getElements();
   double *x = fp_data->x_ip;

   *is_feasible = TRUE;
   /* some int variable is non-integral */
   /* is not possible, since this function is called after rounding */

   /* check feasibility of constraints */
   for (i=0;i<m;i++) {
      Ractivity = 0;
      c=0;			/* column */
      for (j=r_matbeg[i];j<r_matbeg[i]+r_matlen[i];j++) {
         c=r_matind[j];
         Ractivity += x[c]*r_matval[j];
      }
      //      printf("Ractivity[%d] = \t%f\n",i,Ractivity);
      if (Ractivity>r_up[i]+lpetol || Ractivity<r_low[i]-lpetol) {
         /* constraint infeasibility is possible since we call this func. after
            rounding */
         *is_feasible = FALSE;
         //printf("constraint %d activity = %f, down = %g, up = %g\n",
	 //i, Ractivity, r_low[i], r_up[i]);
         break;
      }
   }

   return termcode;
}

/*===========================================================================*/
int fp_initialize_lp_solver(lp_prob *p, LPdata *new_lp_data, FPdata *fp_data)
{
   /*
      create a copy of lp_data into new_lp_data
      for general mixed int programs, we will have to add 2 new vars for each
      non-binary integer var. (x_j+ and x_j-)
   */

   /* first create an exact copy of lp_data */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   LPdata *lp_data  = p->lp_data;
   new_lp_data->lpetol = lp_data->lpetol;
   int n = lp_data->n;
   int m = lp_data->m;
   int i, k, *outrhsind;
   //int *rstat,*cstat;

   double one=1.0;
   char sense='G';
   char where_to_move='E';	/* redundant */
   int col_number = n;
   int *rmatbeg = (int *) malloc(2*ISIZE);
   int *cmatbeg = (int *) malloc(2*ISIZE);
   int *rmatind = (int *) malloc(3*ISIZE);
   double *rmatval = (double *) malloc(3*DSIZE);
   int *cmatind = NULL;
   double *cmatval = NULL;
   double rhs;
   double lb, ub;
   double lpetol = lp_data->lpetol;
   double *lp_lb, *lp_ub, *fp_obj;
   double norm_c = 0;
   double *mip_obj = fp_data->mip_obj;
   int verbosity = fp_data->verbosity;
   int *index_list;
   int fp_max_length_cuts = 1;
   row_data *rows = lp_data->rows;

   /* used because we can not call si directly */
   copy_lp_data(lp_data,new_lp_data);
#ifdef __OSI_CLP__
   new_lp_data->si->setupForRepeatedUse(3,0);
#endif

#ifdef COMPILE_IN_LP
   if(p->mip->matbeg){
     double mat_den = (1.0)*p->mip->nz/(p->mip->m * p->mip->n + 1);
#ifdef __OSI_CLP__
     if(p->mip->nz > 1e5 && mat_den > 0.01){
       new_lp_data->si->setupForRepeatedUse(0,0);
     }
#endif
   }
#endif

   lp_lb = new_lp_data->lb;
   lp_ub = new_lp_data->ub;

   /* delete cuts that are long as they slow down the lp */
   outrhsind = (int *)calloc(m, ISIZE);
   k = 0;

#ifdef COMPILE_IN_LP
   if(p->bc_level < 1 && p->mip->mip_inf && p->mip->mip_inf->cont_var_num <= 0){
      fp_max_length_cuts = 100;
   }
#endif

   for (i = p->base.cutnum; i < m; i++){
      if (((int *)rows[i].cut->coef)[0] > fp_max_length_cuts) {
         outrhsind[k] = i;
         k++;
      }
   }
   PRINT(verbosity, 5, ("feasibility pump: cuts discarded = %d\n", k));
   delete_rows_with_ind(new_lp_data, k, outrhsind);
   m -= k;
   //   printf("m: %i \n",m);
   /* set up fp_data */
   fp_data->alpha           = 0.8;
   fp_data->alpha_decr      = 0.7;
   fp_data->n0 = fp_data->n = n;

   fp_data->m0              = m;
   fp_data->iter            = 0;

   /* count how many binary variables */
   fp_data->fp_vars         = (FPvars **) malloc(sizeof(FPvars *)*n);
   fp_data->x_ip            = (double *) calloc(n,DSIZE);
   fp_data->x_lp            = (double *) calloc(n,DSIZE);
   fp_data->index_list      = (int *)    calloc(n,DSIZE);
   fp_data->x_bar_ind       = (int **)   calloc(p->par.fp_max_cycles,
                                                sizeof(int*));
   fp_data->x_bar_val       = (double **)calloc(p->par.fp_max_cycles,
                                                sizeof(double*));
   fp_data->x_bar_len       = (int *)    calloc(p->par.fp_max_cycles,ISIZE);
   fp_data->alpha_p         = (double *) malloc(p->par.fp_max_cycles*DSIZE);
   FPvars **fp_vars         = fp_data->fp_vars;
   fp_data->numNonBinInts   = 0;
   fp_data->numInts         = 0;

   index_list = fp_data->index_list;
   for (i=0;i<n;i++) {
      index_list[i]=i;
      fp_vars[i] = (FPvars *)malloc(sizeof(FPvars));
      if (lp_data->vars[i]->is_int) {
         fp_data->numInts++;
         fp_vars[i]->is_int = TRUE;
         if (lp_data->vars[i]->lb<-lpetol||lp_data->vars[i]->ub>1+lpetol) {
            fp_vars[i]->is_bin = FALSE;
            fp_data->numNonBinInts++;
         }
         else {
            fp_vars[i]->is_bin = TRUE;
         }
      } else {
         fp_vars[i]->is_int = fp_vars[i]->is_bin = FALSE;
      }
      /* calculate ||C|| */
      norm_c += mip_obj[i]*mip_obj[i];
   }

   norm_c = sqrt(norm_c);
   PRINT(verbosity, 20, ("fp: norm_c = %f\n",norm_c));

   fp_data->n       = n+fp_data->numNonBinInts;
   fp_data->m       = m+2*fp_data->numNonBinInts;
   fp_data->obj     = (double *)malloc(fp_data->n*DSIZE);
   new_lp_data->x   = (double *)calloc(fp_data->n,DSIZE);
   memcpy(fp_data->x_lp,p->lp_data->x,DSIZE*n);

   if (norm_c>lpetol) {
      for (i=0;i<n;i++) {
         mip_obj[i] = mip_obj[i]/norm_c;
      }
   }

   /* load basis */
   //rstat = (int *) malloc(m * ISIZE);
   //cstat = (int *) malloc(n * ISIZE);

   //get_basis(lp_data,cstat,rstat);
   //load_basis (new_lp_data,cstat,rstat);

   //FREE(rstat);
   //FREE(cstat);

   /* add 1 columns and 2 rows for each nonBinary Integer */
   /*
    * min d_i
    * s.t.
    * d_i - x_i >= -x_i^0
    * d_i + x_i >=  x_i^0
    */
   rmatbeg[0] =  0;
   rmatbeg[1] =  2;
   cmatbeg[0] =  0;
   cmatbeg[1] =  0;
   rmatval[0] =  1.0;
   lb         = -SYM_INFINITY;
   ub         =  SYM_INFINITY;
   fp_obj     =  fp_data->obj;

   for (i=0;i<n;i++) {
      if (fp_vars[i]->is_int && !fp_vars[i]->is_bin) {
         /* add d_i */
         add_cols(new_lp_data, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub,
               &where_to_move);
         fp_vars[i]->xplus = col_number;

         /* now add two rows */
         /* d_i - x_i >= -x_i^0 */
         rhs        = -1*lp_data->x[i];
         rmatind[0] =  col_number;
         rmatind[1] =  i;
         rmatval[1] = -1.0;
         add_rows(new_lp_data, 1, 2, &rhs, &sense, rmatbeg, rmatind, rmatval);

         /* d_i - x_i >= -x_i^0 */
         rhs = lp_data->x[i];
         rmatval[1] = 1.0;
         add_rows(new_lp_data, 1, 2, &rhs, &sense, rmatbeg, rmatind, rmatval);

         fp_obj[col_number] = 1.0;
         col_number++;
      }
   }

   /* used by change_rhs */
   new_lp_data->tmp.c = (char *)malloc(2*CSIZE);
   new_lp_data->tmp.d = (double *)malloc(DSIZE*n);
   new_lp_data->tmp.i1 = (int *)malloc(ISIZE*n);

   FREE(rmatval);
   FREE(rmatind);
   FREE(cmatbeg);
   FREE(rmatbeg);
   FREE(outrhsind);

   return termcode;
}

/*===========================================================================*/
int fp_solve_lp(LPdata *lp_data, FPdata *fp_data, char* is_feasible)
{
   /* construct an lp based on x_ip. solve it. store the result in x_lp */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *objcoeff= fp_data->obj;
   int n = fp_data->n;
   //int iterd;
   int termstatus;
   int i;
   double delta_x;
   double norm = 0;
   FPvars **fp_vars = fp_data->fp_vars;
   double *mip_obj  = fp_data->mip_obj;
   int verbosity = fp_data->verbosity;
   int  *index_list = fp_data->index_list;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   double alpha = fp_data->alpha;
   double one_minus_alpha = 1-fp_data->alpha;
   int n0 = fp_data->n0;
   double *lp_data_x = lp_data->x;
   double etol = lp_data->lpetol;

   is_feasible = FALSE;
   memset ((char *)(objcoeff),0,DSIZE*n);
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         if (fp_vars[i]->is_bin) {
            if (x_ip[i] <= 0.0 + etol && x_ip[i] >= 0.0 - etol) {
               objcoeff[i] = 10.0;
	    } else if (x_ip[i] >= 1.0 - etol && x_ip[i] <= 1.0 + etol ) {
	       objcoeff[i] = -10.0;
            }
         } else {
            objcoeff[i] = 0.0;
            objcoeff[fp_vars[i]->xplus] = 1.0;
         }
      } else {
         objcoeff[i] = 0.0;
      }
      /* calculate ||coeff||, norm is not zero because otherwise x_ip is
       * feasible */
   }

   if (fp_data->iter < 1) {
      norm = 0;
      for (i=0; i < n0; i++) {
         norm += objcoeff[i]*objcoeff[i]; /* stays the same every iteration */
      }
      norm = sqrt(norm);
      fp_data->norm = norm;
   } else {
      norm = fp_data->norm;
   }

   //norm = 0;
   PRINT(verbosity, 15, ("fp: norm = %f\n",norm));
   for (i=0;i<n0;i++) {
      objcoeff[i] =
      one_minus_alpha*objcoeff[i]+alpha*mip_obj[i]*norm;
   }
  /*
   for (i=fp_data->n0;i<fp_data->n;i++) {
      objcoeff[i] = (1-alpha)*objcoeff[i];
   }
   alpha = alpha*fp_data->alpha_decr;
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         lp_data->si->setInteger(i);
      }
   }
   */

   change_objcoeff(lp_data, index_list, &index_list[n-1], objcoeff);
   if (fp_data->iter > 0) {
      termstatus = dual_simplex(lp_data, &fp_data->iterd);
   } else {
      termstatus = initial_lp_solve(lp_data, &fp_data->iterd);
   }

   if (termstatus != LP_OPTIMAL) {
      PRINT(verbosity,0,("Feasibility Pump: Unable to solve LP. Pump malfunction.\n"));
      return FUNCTION_TERMINATED_ABNORMALLY;
   }

   get_x(lp_data);

   delta_x = 0;
   memcpy(x_lp,lp_data_x,DSIZE*n0);

   /*
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         delta_x += fabs(x_lp[i]-x_ip[i]);
      }
   }
   PRINT(verbosity, 15, ("fp: delta_x = %f\n",delta_x));
   */

   return termcode;
}

/*===========================================================================*/
int fp_add_obj_row(LPdata *new_lp_data, int n, const double *obj, double rhs)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   char sense = 'L';
   int *rmatbeg, *rmatind;
   double *rmatval;
   int i, count, nz;
   double lpetol = new_lp_data->lpetol;

   // count non zeros
   // we dont trust p->mip->mip_inf->obj_size because it is the size before
   // preprocessing.
   nz = 0;
   for (i=0;i<n;i++) {
      if (fabs(obj[i])>lpetol) {
         nz++;
      }
   }

   rmatbeg = (int *) malloc(2*ISIZE);
   rmatind = (int *) malloc(nz*ISIZE);
   rmatval = (double *) malloc(nz*DSIZE);

   count = 0;
   for (i=0;i<n;i++) {
      if (fabs(obj[i])>lpetol) {
         rmatval[count] = obj[i];
         rmatind[count] = i;
         count++;
      }
   }
   rmatbeg[0] = 0;
   rmatbeg[1] = nz;
   add_rows(new_lp_data, 1, nz, &rhs, &sense, rmatbeg, rmatind, rmatval);
   FREE(rmatbeg);
   FREE(rmatind);
   FREE(rmatval);
   return termcode;
}

/*===========================================================================*/
int fp_round(lp_prob *p, FPdata *fp_data, LPdata *lp_data)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   int i,j, has_changed;
   int n = fp_data->n0;
   double lpetol = lp_data->lpetol;
   int *tind = lp_data->tmp.i1; /* n */
   double *tx = lp_data->tmp.d; /* n */
   int cnt = 0;
   int *index = fp_data->index_list;
   double **x_bar_val_p = fp_data->x_bar_val;
   double *x_bar_val;
   int **x_bar_ind_p = fp_data->x_bar_ind;
   int *x_bar_ind;
   int *x_bar_len = fp_data->x_bar_len;
   double flip_fraction = fp_data->flip_fraction;
   FPvars **vars = fp_data->fp_vars;
   int fp_iter = fp_data->iter;
   double *alpha_p = fp_data->alpha_p;
   int sos_row_filled_cnt = 0;

   if(fp_data->can_check_sos){
      memset(fp_data->sos_row_filled, 0, CSIZE*p->mip->m);
      //memset(fp_data->sos_var_fixed_zero, 0, CSIZE*p->mip->n);
   }

   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
         /* round x_lp[i] and put into x_ip[i] */
         x_ip[i]=floor(x_lp[i]+0.5);
	 /*
	 if(vars[i]->is_bin && fp_data->can_check_sos && x_ip[i] == 1.0 &&
	    p->mip->mip_inf->cols[i].sos_num){
	    if(fp_data->sos_var_fixed_zero[i]) x_ip[i] = 0;
	    else fp_fix_sos_var(p, fp_data, i);
	 }
	 */
	 if(vars[i]->is_bin && fp_data->can_check_sos && x_ip[i] == 1.0 &&
	    p->mip->mip_inf->cols[i].sos_num){
	    if(!(fp_can_sos_var_fix(p, fp_data, i, &sos_row_filled_cnt))){
	       x_ip[i] = 0.0;
	    }
	 }
      }
      else {
         x_ip[i]=x_lp[i];
      }
   }

   // TODO: make it work for '0'
   //       remove randomness
   while (1) {
      cnt = 0;
      for (i = 0; i < n; i++){
         if (vars[i]->is_int && (x_ip[i] > lpetol || x_ip[i] < -lpetol)){
            tind[cnt] = index[i];
            tx[cnt++] = x_ip[i];
         }
      }
      /* order indices and values according to indices */
      qsort_id(tind, tx, cnt);

      /* go through all 'iter' points and check if x_ip already exists */
      has_changed = TRUE;
      for (i=0; i<fp_iter; i++) {
         //printf("alpha = %f, len = %d\n", alpha_p[i], x_bar_len[i]);
         if (x_bar_len[i] == cnt && alpha_p[i] < 0.08) {
            x_bar_val = x_bar_val_p[i];
            x_bar_ind = x_bar_ind_p[i];
            for (j=0; j<cnt; j++) {
               if (tind[j]!=x_bar_ind[j] || fabs(tx[j]-x_bar_val[j])>lpetol) {
                  break;
               }
            }
            if (j==cnt) {
               PRINT(fp_data->verbosity,5,("fp: same as %d\n",i));
               break; //its same
            }
         }
      }

      if (i<fp_iter) {
         /* flip some vars in x_ip */
	 //if(fp_data->can_check_sos){
	 //  memset(fp_data->sos_row_filled, 0, CSIZE*p->mip->m);
	 //  sos_row_filled_cnt = 0;
	 //}

         int num_flipped = 0;

         has_changed = FALSE;
         PRINT(fp_data->verbosity,5,("fp: flipping\n"));

	 for (j=0; j<n; j++) {
	    if (vars[j]->is_bin) {

	       if (CoinDrand48()<flip_fraction) {
		  x_ip[j] = 1-x_ip[j];
		  num_flipped++;
	       }
	       // if(fp_data->can_check_sos && x_ip[j] == 1.0 &&
	       //  p->mip->mip_inf->cols[j].sos_num){
		  //if(!(fp_can_sos_var_fix(p, fp_data, j, &sos_row_filled_cnt))){
		  //   x_ip[j] = 0.0;
		  // }
	       // }
	    } else if (vars[j]->is_int) {
	       if (CoinDrand48()<flip_fraction) {
		  x_ip[j] = floor(x_lp[j]) +
		     floor(ceil(x_lp[j]) - x_lp[j] + 0.5); /*round and flip*/
	       }
	    }
	 }

	 PRINT(fp_data->verbosity,5,("fp: flipping %d\n", num_flipped));
         if (num_flipped==0) {
            // TODO: dont know what to do
            break;
         }
      } else {
         break;
      }
   }

   /*
   int k;
   if(fp_data->can_check_sos && p->mip->mip_inf->binary_sos_row_num > sos_row_filled_cnt){
      int fix_col = 0;
      int row_ind = 0;
      for(k = 0; k < p->mip->m; k++){
	 if(p->mip->mip_inf->rows[k].is_sos_row &&
	    !(fp_data->sos_row_filled[k])){
	    fix_col = p->mip->row_matind[p->mip->row_matbeg[k]];
	    for(j = p->mip->matbeg[fix_col]; j < p->mip->matbeg[fix_col + 1];
		j++){
	       row_ind = p->mip->matind[j];
	       if(p->mip->mip_inf->rows[row_ind].is_sos_row){
		  fp_data->sos_row_filled[row_ind] = TRUE;
		  sos_row_filled_cnt++;
	       }
	    }
	    x_ip[fix_col] = 1.0;
	    if(sos_row_filled_cnt >= p->mip->mip_inf->binary_sos_row_num){
	       break;
	    }
	 }
      }
   }
   */

   if (has_changed==TRUE || fp_data->alpha>0) {
      fp_data->x_bar_ind[fp_iter] = (int *)malloc(ISIZE*cnt);
      fp_data->x_bar_val[fp_iter] = (double *)malloc(DSIZE*cnt);
      x_bar_len[fp_iter] = cnt;
      memcpy(fp_data->x_bar_ind[fp_iter],tind,ISIZE*cnt);
      memcpy(fp_data->x_bar_val[fp_iter],tx,DSIZE*cnt);
      fp_data->alpha = fp_data->alpha*fp_data->alpha_decr;
      if (fp_data->alpha<0.08) {
         fp_data->alpha = 0;
      }
      fp_data->alpha_p[fp_iter] = fp_data->alpha;
   } else {
      x_bar_len[fp_iter] = -1;
   }
   return termcode;
}
/*===========================================================================*/

int fp_fix_sos_var(lp_prob *p, FPdata *fp_data, int ind)
{

   int k, j, row_ind, col_ind;
   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      for(j = p->mip->row_matbeg[row_ind + 1] - 1; j >= p->mip->row_matbeg[row_ind] ; j--){
	 col_ind = p->mip->row_matind[j];
	 if(col_ind <= ind) break;
	 else fp_data->sos_var_fixed_zero[col_ind] = TRUE;
      }
   }

   return 0;
}

/*===========================================================================*/

int fp_can_sos_var_fix(lp_prob *p, FPdata *fp_data, int ind, int *filled_row_cnt)
{
   int k, row_ind;

   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      if(p->mip->mip_inf->rows[row_ind].is_sos_row){
	 if(fp_data->sos_row_filled[row_ind]){
	    return FALSE;
	 }
      }
   }
   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      if(p->mip->mip_inf->rows[row_ind].is_sos_row){
	 fp_data->sos_row_filled[row_ind] = TRUE;
	 (*filled_row_cnt)++;
      }
   }

   return TRUE;
}
/*===========================================================================*/
int fp_should_call_fp(lp_prob *p, int branching, int *should_call,
      char is_last_iter)
{
   int        termcode = FUNCTION_TERMINATED_NORMALLY;

   *should_call = FALSE;
   if (is_last_iter==FALSE || (p->has_ub && p->lp_stat.fp_calls > 100)){
			       //(p->ub-p->lp_data->objval)/(fabs(p->ub)+0.0001)*100 >
			       //2*p->par.fp_min_gap)){
      return termcode;
   }

   int fp_freq_base = p->bc_level;
#ifdef COMPILE_IN_LP
   //   fp_freq_base = p->tm->stat.analyzed - 1;
#endif

   int orig_fp_freq = p->par.fp_frequency;
   if(!p->has_ub && p->lp_stat.fp_calls < 3 &&
      p->lp_stat.lp_total_iter_num/(p->lp_stat.lp_calls -
				    p->lp_stat.str_br_lp_calls -
				    p->lp_stat.fp_lp_calls + 1) > 1000){
      p->par.fp_frequency = 5;
   }

   if (p->par.fp_enabled>0 && !branching) {
      if (p->par.fp_enabled == SYM_FEAS_PUMP_REPEATED &&
	  (fp_freq_base)%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if (p->has_ub==FALSE && p->par.fp_enabled==SYM_FEAS_PUMP_TILL_SOL
            && p->bc_level%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if (  (p->has_ub==FALSE||
		    (p->ub-p->lp_data->objval)/(fabs(p->ub)+0.0001)*100>
		    p->par.fp_min_gap) &&
		   (p->comp_times.fp < p->par.fp_max_initial_time) &&// ||
		   //p->lp_stat.fp_lp_total_iter_num <
		   //0.025*(p->lp_stat.lp_total_iter_num +
		   //   p->lp_stat.fp_lp_total_iter_num)) &&
		   //p->comp_times.fp < 0.025*p->tt) &&
		   //p->comp_times.fp < 0.5*p->tt && //menal, also index->level
		   fp_freq_base%p->par.fp_frequency == 0 ) {
         *should_call = TRUE;
      }
   }

   if(p->bc_level < 1 && p->lp_stat.fp_calls > 0 &&
      p->comp_times.fp >= 0.5*p->par.fp_time_limit){
      *should_call = FALSE;
   }else if (!should_call){
      if(p->bc_level > 0 && !p->has_ub && p->lp_stat.fp_calls <= 3){
	 *should_call = TRUE;
      }
   }

   p->par.fp_frequency = orig_fp_freq;

   if (*should_call == TRUE) {
      p->lp_stat.num_fp_calls_in_path++;
   }
   return termcode;
}

/*===========================================================================*/
// menal - adapted from cbc
// See if rounding will give solution

int round_solution(lp_prob *p, double *solutionValue, double *betterSolution)
{

   int numberColumns = p->mip->n;
   int numberRows = p->mip->m;
   int nz = p->mip->nz;
   int returnCode = 0, numberIntegers = 0;
   double primalTolerance = p->lp_data->lpetol,
      integerTolerance = primalTolerance;
   double *lower, *upper, *solution, *objective;
   double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
   double newSolutionValue = direction*p->lp_data->objval;
   double *element, *elementByRow;
   int * integerVariable, row_ind, elem_ind;
   int *row, *column, *columnStart, *rowStart, *columnLength, *rowLength;
   int i, j;

   if(!(p->mip->matbeg)){
     return returnCode;
   }

   get_bounds(p->lp_data);
   get_x(p->lp_data);

   lower = p->lp_data->lb;
   upper = p->lp_data->ub;
   solution = p->lp_data->x;

   element = p->mip->matval;
   row = p->mip->matind;
   columnStart = p->mip->matbeg;
   objective = p->mip->obj;

   columnLength = p->mip->col_lengths;

   if(!columnLength){
      columnLength=(p->mip->col_lengths = (int *)calloc(numberColumns,ISIZE));
      elementByRow = (p->mip->row_matval = (double *)malloc(nz*DSIZE));
      column = (p->mip->row_matind = (int *)malloc(nz*ISIZE));
      rowStart = (p->mip->row_matbeg = (int *)malloc((numberRows+1)*ISIZE));
      rowLength = (p->mip->row_lengths = (int *)calloc(numberRows,ISIZE));

      /* first get row legths */
      for(i = 0; i < numberColumns; i++){
	 /* get orig indices here */
	 for(j = columnStart[i]; j < columnStart[i+1]; j++){
	    rowLength[row[j]]++;
	 }
	 columnLength[i] = columnStart[i+1] - columnStart[i];
      }

      rowStart[0] = 0;

      /* fill in matbegs */
      for(i = 0; i < numberRows; i++){
	 rowStart[i + 1] = rowStart[i] + rowLength[i];
      }

      /* get matrix, change 'G' rows to 'L'*/
      for(i = 0; i < numberColumns; i++){
	 for(j = columnStart[i]; j < columnStart[i+1]; j++){
	    row_ind = row[j];
	    elem_ind = rowStart[row_ind];
	    column[elem_ind] = i;

	    elementByRow[elem_ind] = element[j];
	    rowStart[row_ind] = elem_ind + 1;
	 }
      }

      for(i = 0; i < numberRows; i++){
	 rowStart[i] -= rowLength[i];
      }
   }else{

      elementByRow = p->mip->row_matval;
      column = p->mip->row_matind;
      rowStart = p->mip->row_matbeg;
      rowLength = p->mip->row_lengths;
   }


   const double * rowUpper = p->lp_data->si->getRowUpper();
   const double * rowLower = p->lp_data->si->getRowLower();

   integerVariable = new int[numberColumns];

   for (i = 0; i<numberColumns; i++){
      if (p->mip->is_int[i]){
	 integerVariable[numberIntegers++] = i;
      }
   }

   // Get solution array for heuristic solution

   double * newSolution = new double [numberColumns];
   memcpy(newSolution,solution,numberColumns*sizeof(double));

   double * rowActivity = new double[numberRows];
   memset(rowActivity,0,numberRows*sizeof(double));
   for (i=0;i<numberColumns;i++) {
      int j;
      double value = newSolution[i];
      if (value) {
	 for (j=columnStart[i];
	      j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    //	printf("rowind %i: %i \n", j, iRow);
	    //	if(j < 5){
	    //	printf("element %i: %f \n", j, element[j]);
	    //	}
	    rowActivity[iRow] += value*element[j];
	 }
      }
   }
   // check was feasible - if not adjust (cleaning may move)
   for (i=0;i<numberRows;i++) {
      if(rowActivity[i]<rowLower[i]) {
	 //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
	 rowActivity[i]=rowLower[i];
      } else if(rowActivity[i]>rowUpper[i]) {
	 //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
	 rowActivity[i]=rowUpper[i];
      }
   }
   for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	 double below = floor(value);
	 double newValue=newSolution[iColumn];
	 double cost = direction * objective[iColumn];
	 double move;
	 if (cost>0.0) {
	    // try up
	    move = 1.0 -(value-below);
	 } else if (cost<0.0) {
	    // try down
	    move = below-value;
	 } else {
	    // won't be able to move unless we can grab another variable
	    // just for now go down
	    move = below-value;
	 }
	 newValue += move;
	 newSolution[iColumn] = newValue;
	 newSolutionValue += move*cost;
	 int j;
	 for (j=columnStart[iColumn];
	      j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    rowActivity[iRow] += move*element[j];
	 }
      }
   }

   double penalty=0.0;

   // see if feasible
   for (i=0;i<numberRows;i++) {
      double value = rowActivity[i];
      double thisInfeasibility=0.0;
      if (value<rowLower[i]-primalTolerance)
	 thisInfeasibility = value-rowLower[i];
      else if (value>rowUpper[i]+primalTolerance)
	 thisInfeasibility = value-rowUpper[i];
      if (thisInfeasibility) {
	 // See if there are any slacks I can use to fix up
	 // maybe put in coding for multiple slacks?
	 double bestCost = 1.0e50;
	 int k;
	 int iBest=-1;
	 double addCost=0.0;
	 double newValue=0.0;
	 double changeRowActivity=0.0;
	 double absInfeasibility = fabs(thisInfeasibility);
	 for (k=rowStart[i];k<rowStart[i]+rowLength[i];k++) {
	    int iColumn = column[k];
	    if (columnLength[iColumn]==1) {
	       double currentValue = newSolution[iColumn];
	       double elementValue = elementByRow[k];
	       double lowerValue = lower[iColumn];
	       double upperValue = upper[iColumn];
	       double gap = rowUpper[i]-rowLower[i];
	       double absElement=fabs(elementValue);
	       if (thisInfeasibility*elementValue>0.0) {
		  // we want to reduce
		  if ((currentValue-lowerValue)*absElement>=absInfeasibility) {
		     // possible - check if integer
		     double distance = absInfeasibility/absElement;
		     double thisCost = -direction*objective[iColumn]*distance;
		     if (p->mip->is_int[iColumn]) {
			distance = ceil(distance-primalTolerance);
			if (currentValue-distance>=lowerValue-primalTolerance) {
			   if (absInfeasibility-distance*absElement< -gap-primalTolerance)
			      thisCost=1.0e100; // no good
			   else
			      thisCost = -direction*objective[iColumn]*distance;
			} else {
			   thisCost=1.0e100; // no good
			}
		     }
		     if (thisCost<bestCost) {
			bestCost=thisCost;
			iBest=iColumn;
			addCost = thisCost;
			newValue = currentValue-distance;
			changeRowActivity = -distance*elementValue;
		     }
		  }
	       } else {
		  // we want to increase
		  if ((upperValue-currentValue)*absElement>=absInfeasibility) {
		     // possible - check if integer
		     double distance = absInfeasibility/absElement;
		     double thisCost = direction*objective[iColumn]*distance;
		     if (p->mip->is_int[iColumn]) {
			distance = ceil(distance-primalTolerance);
		//assert (currentValue-distance<=upperValue+primalTolerance);
			if (absInfeasibility-distance*absElement< -gap-primalTolerance)
			   thisCost=1.0e100; // no good
			else
			   thisCost = direction*objective[iColumn]*distance;
		     }
		     if (thisCost<bestCost) {
			bestCost=thisCost;
			iBest=iColumn;
			addCost = thisCost;
			newValue = currentValue+distance;
			changeRowActivity = distance*elementValue;
		     }
		  }
	       }
	    }
	 }
	 if (iBest>=0) {
	    /*printf("Infeasibility of %g on row %d cost %g\n",
	      thisInfeasibility,i,addCost);*/
	    newSolution[iBest]=newValue;
	    thisInfeasibility=0.0;
	    newSolutionValue += addCost;
	    rowActivity[i] += changeRowActivity;
	 }
	 penalty += fabs(thisInfeasibility);
      }
   }

   // Could also set SOS (using random) and repeat
   if (!penalty) {
      // See if we can do better
      //seed_++;
      //CoinSeedRandom(seed_);
      // Random number between 0 and 1.
      double randomNumber = CoinDrand48();
      int iPass;
      int start[2];
      int end[2];
      int iRandom = (int) (randomNumber*((double) numberIntegers));
      start[0]=iRandom;
      end[0]=numberIntegers;
      start[1]=0;
      end[1]=iRandom;
      for (iPass=0;iPass<2;iPass++) {
	 int i;
	 for (i=start[iPass];i<end[iPass];i++) {
	    int iColumn = integerVariable[i];
	    //double value=newSolution[iColumn];
	    //assert (fabs(floor(value+0.5)-value)<integerTolerance);
	    double cost = direction * objective[iColumn];
	    double move=0.0;
	    if (cost>0.0)
	       move = -1.0;
	    else if (cost<0.0)
	       move=1.0;
	    while (move) {
	       bool good=true;
	       double newValue=newSolution[iColumn]+move;
	       if (newValue<lower[iColumn]-primalTolerance||
		   newValue>upper[iColumn]+primalTolerance) {
		  move=0.0;
	       } else {
		  // see if we can move
		  int j;
		  for (j=columnStart[iColumn];
		       j<columnStart[iColumn]+columnLength[iColumn];j++) {
		     int iRow = row[j];
		     double newActivity = rowActivity[iRow] + move*element[j];
		     if (newActivity<rowLower[iRow]-primalTolerance||
			 newActivity>rowUpper[iRow]+primalTolerance) {
			good=false;
			break;
		     }
		  }
		  if (good) {
		     newSolution[iColumn] = newValue;
		     newSolutionValue += move*cost;
		     int j;
		     for (j=columnStart[iColumn];
			  j<columnStart[iColumn]+columnLength[iColumn];j++) {
			int iRow = row[j];
			rowActivity[iRow] += move*element[j];
		     }
		  } else {
		     move=0.0;
		  }
	       }
	    }
	 }
      }
      if (newSolutionValue < *solutionValue) {
	 // paranoid check
	 memset(rowActivity,0,numberRows*sizeof(double));
	 for (i=0;i<numberColumns;i++) {
	    int j;
	    double value = newSolution[i];
	    if (value) {
	       for (j=columnStart[i];
		    j<columnStart[i]+columnLength[i];j++) {
		  int iRow=row[j];
		  rowActivity[iRow] += value*element[j];
	       }
	    }
	 }
	 // check was approximately feasible
	 bool feasible=true;
	 for (i=0;i<numberRows;i++) {
	    if(rowActivity[i]<rowLower[i]) {
	       if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
		  feasible = false;
	    } else if(rowActivity[i]>rowUpper[i]) {
	       if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
		  feasible = false;
	    }
	 }
	 if (feasible) {
	    // new solution
	    memcpy(betterSolution, newSolution, numberColumns*DSIZE);
	    *solutionValue = newSolutionValue;
	    //printf("** Solution of %g found by rounding\n",newSolutionValue);
	    returnCode=1;
	 } else {
	    // Can easily happen
	    //printf("Debug CbcRounding giving bad solution\n");
	 }
      }
   }
  delete [] integerVariable;

  delete [] newSolution;
  delete [] rowActivity;
  return returnCode;
}

/*===========================================================================*/
/* --menal
  adapted from cbc-disabled
*/
/*===========================================================================*/
int local_search(lp_prob *p, double *solutionValue, double *colSolution,
		 double *betterSolution)
{

   LPdata *lp_data = p->lp_data;
   int numberColumns = p->mip->n;
   int numberRows = p->mip->m, nz = p->mip->nz;
   int returnCode = 0, numberIntegers = 0;
   double primalTolerance = lp_data->lpetol;
   double *solution = colSolution, *objective;
   double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
   double newSolutionValue = direction*(*solutionValue);
   double *element, *elementByRow;
   int * integerVariable, *rowStart;
   int *row, *columnStart, *columnLength, *column, *rowLength;
   int i, j, row_ind, elem_ind;


  element = p->mip->matval;
  row = p->mip->matind;
  columnStart = p->mip->matbeg;
  objective = p->mip->obj;

  columnLength = p->mip->col_lengths;

  if(!columnLength){
     columnLength=(p->mip->col_lengths = (int *)calloc(numberColumns,ISIZE));
     elementByRow = (p->mip->row_matval = (double *)malloc(nz*DSIZE));
     column = (p->mip->row_matind = (int *)malloc(nz*ISIZE));
     rowStart = (p->mip->row_matbeg = (int *)malloc((numberRows+1)*ISIZE));
     rowLength = (p->mip->row_lengths = (int *)calloc(numberRows,ISIZE));

     /* first get row legths */
     for(i = 0; i < numberColumns; i++){
	/* get orig indices here */
	for(j = columnStart[i]; j < columnStart[i+1]; j++){
	   rowLength[row[j]]++;
	}
	columnLength[i] = columnStart[i+1] - columnStart[i];
     }

     rowStart[0] = 0;

     /* fill in matbegs */
     for(i = 0; i < numberRows; i++){
	rowStart[i + 1] = rowStart[i] + rowLength[i];
     }

     /* get matrix, change 'G' rows to 'L'*/
     for(i = 0; i < numberColumns; i++){
	for(j = columnStart[i]; j < columnStart[i+1]; j++){
	   row_ind = row[j];
	   elem_ind = rowStart[row_ind];
	   column[elem_ind] = i;

	   elementByRow[elem_ind] = element[j];
	   rowStart[row_ind] = elem_ind + 1;
	}
     }

     for(i = 0; i < numberRows; i++){
	rowStart[i] -= rowLength[i];
     }
  }else{

     elementByRow = p->mip->row_matval;
     column = p->mip->row_matind;
     rowStart = p->mip->row_matbeg;
     rowLength = p->mip->row_lengths;
  }


  const double * rowUpper = p->lp_data->si->getRowUpper();
  const double * rowLower = p->lp_data->si->getRowLower();

   integerVariable = new int[numberColumns];

   for (i = 0; i<numberColumns; i++){
      if (p->mip->is_int[i]){
	 integerVariable[numberIntegers++] = i;
      }
   }

  // Column copy
  /*
  const double * element = matrix.getElements();
  const int * row = matrix.getIndices();
  const CoinBigIndex * columnStart = matrix.getVectorStarts();
  const int * columnLength = matrix.getVectorLengths();
  */

  // Get solution array for heuristic solution
  double * newSolution = new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  // way is 1 if down possible, 2 if up possible, 3 if both possible
  char * way = new char[numberIntegers];
  // corrected costs
  double * cost = new double[numberIntegers];
  // for array to mark infeasible rows after iColumn branch
  char * mark = new char[numberRows];
  memset(mark,0,numberRows);
  // space to save values so we don't introduce rounding errors
  double * save = new double[numberRows];

  // clean solution
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];

    // get original bounds
    //    double originalLower = lp_data->vars[iColumn]->lb; //p->mip->lb[iColumn];
    // double originalUpper = lp_data->vars[iColumn]->ub; //p->mip->ub[iColumn];

    double originalLower = p->mip->lb[iColumn];
    double originalUpper = p->mip->ub[iColumn];

    double value=newSolution[iColumn];

    if (value<originalLower) {
       value=originalLower;
       newSolution[iColumn]=value;
    } else if (value>originalUpper) {
       value=originalUpper;
       newSolution[iColumn]=value;
    }

    double nearest=floor(value+0.5);
    //assert(fabs(value-nearest)<10.0*primalTolerance);
    value=nearest;
    newSolution[iColumn]=nearest;
    // if away from lower bound mark that fact
    if (nearest>originalLower) {
      //      used_[iColumn]=1;
    }
    cost[i] = direction*objective[iColumn];
    int iway=0;

    if (value>originalLower+0.5)
      iway = 1;
    if (value<originalUpper-0.5)
       iway |= 2;
    way[i]=(char)iway;
  }
  // get row activities
  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));

  for (i=0;i<numberColumns;i++) {
    int j;
    double value = newSolution[i];
    if (value) {
      for (j=columnStart[i];
	   j<columnStart[i]+columnLength[i];j++) {
	int iRow=row[j];
	rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  // if very infeasible then give up
  bool tryHeuristic=true;
  for (i=0;i<numberRows;i++) {
     if(rowActivity[i]<rowLower[i]) {
      if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      if (rowActivity[i]<rowUpper[i]+10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowUpper[i];
    }
  }
  if (tryHeuristic) {

    // best change in objective
    double bestChange=0.0;

    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];

      double objectiveCoefficient = cost[i];
      int k;
      int j;
      int goodK=-1;
      int wayK=-1,wayI=-1;
      if ((way[i]&1)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] -= element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try down
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (-objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=-1;
		bestChange = -objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (-objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=-1;
		bestChange = -objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if ((way[i]&2)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] += element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try up
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=1;
		bestChange = objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=1;
		bestChange = objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if (goodK>=0) {
	// we found something - update solution
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayI * element[j];
	}
	newSolution[iColumn] += wayI;
	int kColumn = integerVariable[goodK];
	for (j=columnStart[kColumn];
	     j<columnStart[kColumn]+columnLength[kColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayK * element[j];
	}
	newSolution[kColumn] += wayK;
	// See if k can go further ?
	// get original bounds
	double originalLower = p->mip->lb[kColumn];
	double originalUpper = p->mip->ub[kColumn];

	double value=newSolution[kColumn];
	int iway=0;
	if (value>originalLower+0.5)
	  iway = 1;
	if (value<originalUpper-0.5)
	  iway |= 2;
	way[goodK]=(char)iway;
      }
    }
    if (bestChange+newSolutionValue<*solutionValue) {
       // paranoid check
      memset(rowActivity,0,numberRows*sizeof(double));

      for (i=0;i<numberColumns;i++) {
	int j;
	double value = newSolution[i];
	if (value) {
	  for (j=columnStart[i];
	       j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    rowActivity[iRow] += value*element[j];
	  }
	}
      }
      int numberBad=0;
      double sumBad=0.0;

      // check was approximately feasible
      for (i=0;i<numberRows;i++) {
	 if(rowActivity[i]<rowLower[i]) {
	    sumBad += rowLower[i]-rowActivity[i];
	    if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	       numberBad++;
	 } else if(rowActivity[i]>rowUpper[i]) {
	    sumBad += rowUpper[i]-rowActivity[i];
	    if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
	       numberBad++;
	 }
      }
      if (!numberBad) {
	 for (i=0;i<numberIntegers;i++) {
	    int iColumn = integerVariable[i];
	    // get original bounds
	    double originalLower = p->mip->lb[iColumn];
	    //double originalUpper = p->mip->ub[iColumn];

	    double value=newSolution[iColumn];
	    // if away from lower bound mark that fact
	    if (value>originalLower) {
	       //used_[iColumn]=1;
	    }
	 }
	 // new solution
	 memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	 returnCode=1;
	 *solutionValue = newSolutionValue + bestChange;
      } else {
	 // bad solution - should not happen so debug if see message
	 printf("Local search got bad solution with %d infeasibilities"
		"summing to %g\n",
		numberBad,sumBad);
      }
    }
  }

  delete [] integerVariable;

  delete [] newSolution;
  delete [] rowActivity;
  delete [] way;
  delete [] cost;
  delete [] save;
  delete [] mark;

  return returnCode;
}

/*===========================================================================*/
/*===========================================================================*/












