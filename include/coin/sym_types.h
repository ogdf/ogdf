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

#ifndef _BB_TYPES_H
#define _BB_TYPES_H

#define MAX_CHILDREN_NUM 4
#define MAX_CHANGE_NUM 6

/*===========================================================================*\
 * This file contains those type definitions which are used in more than one
 * process of the black box model.
\*===========================================================================*/

#if !defined (_MSC_VER)
#include <unistd.h>            /* this defines sleep() */
#endif

#include "sym_proto.h"

typedef struct LP_SOL{
   int            lp;          /* the tid of the lp process asssociated with
				  the current solution */
   int            has_sol;     /* indicates whether a feasible solution
				  found*/
   int            xlength;     /* the number of nonzeros in the lp solution
				  currently being processed*/
   int            xlevel;      /* the level at which the current solution was
				  generated*/
   int            xindex;
   int            xiter_num;
   int            max_sol_length;
   int           *xind;        /* the indices of the nonzeros in the current
				  solution*/
   double        *xval;        /* the values of the nonzeros in the current
				  solution*/
   double         objval;      /* the objective function value of the current
				  relaxation*/
   double         lpetol;
}lp_sol;

typedef struct BASE_DESC{
   int            varnum;
   int           *userind;
#if 0
   double        *lb;          /* even if there are global lb and ub, we */
   double        *ub;          /* fill these arrays out */
#endif
   int            cutnum;
}base_desc;

/*===========================================================================*\
 * Stores the data corresponding to a particular cut
\*===========================================================================*/

typedef struct CUT_DATA{
   int            size;        /* the size of the coef array */
   char          *coef;        /* an array which contains the data necessary to
				  construct the cut -- it is stored in a
				  packed form. The types of the cut tells how
				  to "unpack" it */
   double         rhs;         /* the right hand side for the constraint*/
   double         range;
   char           type;        /* the type of cut */
   char           sense;       /* sense of the cut constraint */
   char           deletable;   /* whether or not this cut should be removed
				  from the LP after being added */
   int            branch;      /* shows whether we can branch on its cut if
				  this row becomes slack */
   int            name;        /* internal to the BB. The identifier of the
				  cut. >=0 if exists, -1 if does not exist yet,
				  but the cuT is sent to the cutpool, -2 if
				  no name & no pool */
}cut_data;

typedef struct ROW_DATA{
   cut_data      *cut;
   int            ineff_cnt;
   int            eff_cnt;
   char           free;
   char           deletable;
}row_data;

typedef struct WAITING_ROW{
   int            source_pid;
   cut_data      *cut;
   int           *matind;
   double        *matval;
   int            nzcnt;
   double         violation;
}waiting_row;

/*===========================================================================*\
 * The following three definitions are used to describe the search tree
 * nodes.
\*===========================================================================*/

typedef struct VAR_DESC{
   int            userind;
   int            colind;
   double         lb;     /* lb on var before we start processing lp_prob */
   double         ub;     /* ub on var before we start processing lp_prob */
   double         new_lb; /* lb may change due to rc fixing or cuts */
   double         new_ub; /* ub may change due to rc fixing or cuts */
   char           is_int; /* whether or not the variable is integer */
}var_desc;

/*========================================================================*\
 * changes in bounds of variables at a node are stored here.
 * Only the changes pertaining to this node are stored here. changes in
 * parents are not stored with children. Variable bounds can changed due to
 * reduced cost fixing, cuts etc.
 \*========================================================================*/
typedef struct BOUNDS_CHANGE_DESC{
   int                 num_changes; /* how many bounds changed */
   int                *index;       /* max size 2*n */
   char               *lbub;        /* ub or lb? */
   double             *value;       /* new bound value */
}bounds_change_desc;


typedef struct BRANCH_DESC{
   int            name;        /* the userind/cut name depending on the type */
   char           type;        /* description of the child. All of them are
				  natural for branching cuts. For branching
				  variables they should be interpreted as if
				  we were adding a cut with a single variable
				  on the left hand side */
   char           sense;
   double         rhs;
   double         range;
   int            branch;
}branch_desc;

typedef struct ARRAY_DESC{
   char           type;        /* NO_DATA_STORED, EXPLICIT_LIST, WRT_PARENT */
   int            size;
   int            added;
   int           *list;
}array_desc;

typedef struct DOUBLE_ARRAY_DESC{
   char           type;        /* NO_DATA_STORED, EXPLICIT_LIST, WRT_PARENT */
   int            size;        /* the size of list, stat */
   int           *list;
   int           *stat;
}double_array_desc;

typedef struct BASIS_DESC{
   char                basis_exists;

   /*========================================================================*\    * Notes:
    *  1) for base...:
    *     if list is non-NULL then it refers to col/row inds, not userinds
    *	    or cut names.
    *  2) the stat field of extra... ponts into the stat field of base...
    *  3) if extra... is EXPLICIT_LIST then the node_desc structure's
    *        cutind/uind fields should be used.
    *  4) EXPLICIT_LIST in uind implies that extravars is explicit,
    *     EXPLICIT_LIST in cutind implies that extrarows is explicit.
   \*========================================================================*/
   double_array_desc   basevars;
   double_array_desc   extravars;
   double_array_desc   baserows;
   double_array_desc   extrarows;
}basis_desc;

typedef struct NODE_DESC{
   /*========================================================================*\
    * The userindices of variables in this node (but not for the base
    * variables); The basis header for this node; The not-yet-permanently-fixed
    * variables (again, no base variable is listed here); and the cuts at
    * this node
   \*========================================================================*/
   array_desc     uind;
   basis_desc     basis;
   array_desc     not_fixed;
   int            nf_status;   /* NF_CHECK_ALL, NF_CHECK_AFTER_LAST,
				  NF_CHECK_UNTIL_LAST, NF_CHECK_NOTHING */
   array_desc     cutind;
#if defined(COMPILING_FOR_LP) || defined(COMPILING_FOR_MASTER) || defined(COMPILE_IN_LP)
   cut_data     **cuts;        /* this is not used in TM anyway. */
#endif

   bounds_change_desc *bnd_change; /* changes in variable bounds that happen
                                     during the processing of the node */

   /* Any additional info the user might want to pass */
   int           desc_size;
   char         *desc;
}node_desc;

typedef struct BRANCH_OBJ{
   char          type;         /* Type of the candidate */
#if defined(COMPILING_FOR_LP) || defined(COMPILE_IN_LP)
   int           position;     /* The position of the candidate */
   waiting_row  *row;          /* Description of the left hand side; makes
				  sense only for branching cuts */
#endif
   int           child_num;    /* Number of kids */
#if defined(COMPILING_FOR_TM) || defined(COMPILING_FOR_MASTER) || defined(COMPILE_IN_LP)
   int           name;         /* userind for VAR, the index for CUT */
#endif
   double        value;        /* for evaluating pcost */

   /*========================================================================*\
    * Description of the children.
    * All of them are natural for branching cuts.
    * For branching variables they should be interpreted as if we were adding
    * a cut with a single variable on the left hand side
   \*========================================================================*/

#ifdef MAX_CHILDREN_NUM
   char          sense[MAX_CHILDREN_NUM];
   double        rhs[MAX_CHILDREN_NUM];
   double        range[MAX_CHILDREN_NUM];
   int           branch[MAX_CHILDREN_NUM];
#ifdef COMPILE_FRAC_BRANCHING
   int           frac_num[MAX_CHILDREN_NUM];
   int          *frac_ind[MAX_CHILDREN_NUM];
   double       *frac_val[MAX_CHILDREN_NUM];
#endif
#else
   char         *sense;
   double       *rhs;
   double       *range;
   int          *branch;
#ifdef COMPILE_FRAC_BRANCHING
   int          *frac_num;
   int         **frac_ind;
   double      **frac_val;
#endif
#endif

#if defined(COMPILING_FOR_LP) || defined(COMPILE_IN_LP)
   double        lhs;          /* purely for the user */

#ifdef MAX_CHILDREN_NUM
   double        objval[MAX_CHILDREN_NUM];   /* arrays of size 'number' */
   int           termcode[MAX_CHILDREN_NUM];
   int           iterd[MAX_CHILDREN_NUM];
   int           feasible[MAX_CHILDREN_NUM];
   int           is_est[MAX_CHILDREN_NUM];

#else
   double       *objval;   /* arrays of size 'number' */
   int          *termcode;
   int          *iterd;
   int          *feasible;

#endif

#endif
   int          *sol_sizes;
   int         **sol_inds;
   double      **solutions;
#ifdef SENSITIVITY_ANALYSIS
   double      **duals;
#endif

}branch_obj;

/*===========================================================================*/

typedef struct STR_INT{
#ifdef _OPENMP
   char      *str;
#else
   char       str[MAX_LINE_LENGTH +1];
#endif
   int        code;
}str_int;

/*===========================================================================*\
 * This is the time measurement structure for an LP node
\*===========================================================================*/

typedef struct NODE_TIMES{
   double        communication;
   double        lp;
   double        lp_setup;
   double        separation;
   double        fixing;
   double        pricing;
   double        strong_branching;
   double        wall_clock_lp;
   double        ramp_up_tm;
   double        ramp_up_lp;
   double        ramp_down_time;
   double        idle_diving;
   double        idle_node;
   double        idle_names;
   double        idle_cuts;
   double        start_node;
   double        cut_pool;

   /* cuts */
   double        cuts;
   double        gomory_cuts;
   double        knapsack_cuts;
   double        oddhole_cuts;
   double        clique_cuts;
   double        probing_cuts;
   double        mir_cuts;
   double        twomir_cuts;
   double        flowcover_cuts;
   double        rounding_cuts;
   double        lift_and_project_cuts;
   double        landp_cuts;
   double        redsplit_cuts;
   double        dupes_and_bad_coeffs_in_cuts;

   double        fp;                            /* feasibility pump */
   double        primal_heur;                   /* all primal heuristics */
}node_times;

/*===========================================================================*\
 * Here we keep track of the computation time for each of the various
 * parts of the computation
\*===========================================================================*/

typedef struct PROB_TIMES{
   double     readtime;    /* time spent reading in the problem*/
   node_times bc_time;
   double     ub_overhead; /* overhead time used doing the upper bounding */
   double     ub_heurtime; /* actual comp time doing the upper bounding */
   double     lb_overhead; /* overhead time doing the lower bounding */
   double     lb_heurtime; /* actual comp time doing the lower bounding */
}prob_times;

/*===========================================================================*\
 * The bc_node data structure stores the information needed to
 * process a node in the branch and cut tree
\*===========================================================================*/

typedef struct BC_NODE{
   int        bc_index;     /* the identifier of the node */
   int        bc_level;     /* the level in the tree of the node */

   int        lp;           /* the tid of the lp processing the node */
   int        cg;           /* the tid of the cut generator serving the node */
   int        cp;           /* the tid of the cut pool assigned to the node */
   double     lower_bound;  /* the current best objective function value
			       obtained in the subproblem */
   int        update_pc;    /* whether the pseudo cost should be updated after
                               solving the LP */
   double     opt_estimate; /* an estimate of the value of the best feasible
			       solution that could be obtained in this node */
   struct BC_NODE  *parent;
   struct BC_NODE **children;
   branch_obj       bobj;

   node_desc  desc;          /* the description of the node,
			       defined in "sym_types.h" */
   char       node_status;

   int        feasibility_status;
   int        sol_size;
   int       *sol_ind;
   double    *sol;
#ifdef SENSITIVITY_ANALYSIS
   double    *duals;
   double     C_LP;
   double     B_IP;
#endif

#ifdef TRACE_PATH
   char       optimal_path;
#endif

   /* usage of different tools in process chain: fp, cuts, strong branching */
   int         num_cut_iters_in_path;
   int         num_cuts_added_in_path;
   int         num_cuts_slacked_out_in_path;
   double      avg_cuts_obj_impr_in_path;
   double      start_objval;
   double      end_objval;
   char        cuts_tried;
   int         num_str_br_cands_in_path;
   double      avg_br_obj_impr_in_path;
   char        used_str;
   int         t_cnt;

   int         num_fp_calls_in_path;
}bc_node;

/*===========================================================================*\
 * Keeps problem statistics
\*===========================================================================*/

typedef struct PROBLEM_STAT{
   double      root_lb;
   int         cuts_in_pool;
   int         max_depth;          /* keeps track of the deepest level reached
				      in the tree so far */
   int         chains;             /* the number of diving chains */
   int         diving_halts;       /* how many times was an already started
				      dive stopped */
   int         tree_size;          /* number of search tree nodes */
   int         created;            /* the number of created nodes (not
				      necessarily the same as tree_size
				      (trimming...) */
   int         analyzed;           /* the number of analyzed (i.e., CG-LP
				      iteration) nodes (not necessarily same
				      as created, leaves can be cut off
				      without analyzing; trimming) */
   int         leaves_before_trimming;
   int         leaves_after_trimming;
   int         vars_not_priced;    /* How many variables did not price out
				      after the first phase */
   int         nf_status;          /* nf_status of the root node after
				      repricing */
   double      max_vsize;
}problem_stat;

/*===========================================================================*/

typedef struct LP_STAT{
   /* LP solver */
   int         lp_calls;
   int         lp_sols;
   int         lp_total_iter_num; /* number of total simplex iterations */
   int         lp_max_iter_num; /* max of lps' simplex iterations */
   int         str_br_lp_calls; /* no of calls from strong branching */
   int         str_br_bnd_changes; /* no of bounds changed due to strong br */
   int         str_br_nodes_pruned; /* no of nodes pruned by strong br */
   int         str_br_total_iter_num; /* number of total simplex iterations by
					 strong br*/
   int         rel_br_full_solve_num;
   int         rel_br_pc_up_num;
   int         rel_br_up_update;
   int         rel_br_pc_down_num;
   int         rel_br_down_update;
   int         rel_br_impr_num;
   /* cuts */
   int         cuts_generated;
   int         gomory_cuts;
   int         knapsack_cuts;
   int         oddhole_cuts;
   int         clique_cuts;
   int         probing_cuts;
   int         mir_cuts;
   int         twomir_cuts;
   int         flowcover_cuts;
   int         rounding_cuts;
   int         lift_and_project_cuts;
   int         landp_cuts;
   int         redsplit_cuts;

   int         cuts_root;
   int         gomory_cuts_root;
   int         knapsack_cuts_root;
   int         oddhole_cuts_root;
   int         clique_cuts_root;
   int         probing_cuts_root;
   int         mir_cuts_root;
   int         twomir_cuts_root;
   int         flowcover_cuts_root;
   int         rounding_cuts_root;
   int         lift_and_project_cuts_root;
   int         landp_cuts_root;
   int         redsplit_cuts_root;

   int         num_poor_cuts;
   int         num_duplicate_cuts;
   int         num_unviolated_cuts;
   int         cuts_added_to_lps;
   int         cuts_deleted_from_lps;

   int         gomory_calls;
   int         knapsack_calls;
   int         oddhole_calls;
   int         clique_calls;
   int         probing_calls;
   int         mir_calls;
   int         twomir_calls;
   int         flowcover_calls;
   int         rounding_calls;
   int         lift_and_project_calls;
   int         landp_calls;
   int         redsplit_calls;

   /* feasibility pump */
   int         fp_calls;
   int         fp_lp_calls;
   int         fp_num_sols;
   int         fp_poor_sols;
   int         fp_lp_total_iter_num;

   /* usage of different tools in process chain: fp, cuts, strong branching */
   int         num_cut_iters_in_path;
   int         num_cuts_added_in_path;
   int         num_cuts_slacked_out_in_path;
   double      avg_cuts_obj_impr_in_path;
   double      start_objval;
   double      end_objval;
   int         num_str_br_cands_in_path;
   double      avg_br_obj_impr_in_path;

   int         num_fp_calls_in_path;
}lp_stat_desc;


typedef struct RC_DESC{
   int         size;
   int         num_rcs;
   int       **indices;
   double    **values;
   double    **ub;
   double    **lb;
   double     *obj;
   int        *cnt;
}rc_desc;

/*===========================================================================*/
/* Implications */
/*===========================================================================*/
/*===========================================================================*/
typedef struct COL_IMP{

   int col_ind;
  struct COL_IMP *c_next;

}col_imp;

typedef struct IMPVAR{

   int  type; /* ROW, COL */
   int  ind;
   int  fix_type; /*'U', 'L, 'F'
		    for column: improve upper bound, lower bound or fix it
		    for row: all other variables need to be fixed to their
		    'U or 'L' or the row is infea'S'ible
		    however, right now it is same with fix_bounds */
   double val; /* if it is a column impl*/
   struct IMPVAR *right;
   struct IMPVAR *left;

}IMPvar;

typedef struct IMPLIST{

   int size;
   IMPvar * head;
   IMPvar * tail;
}IMPlist;

/*===========================================================================*/
/* Data structure to keep relevant info of a column */
/*===========================================================================*/
typedef struct COLINFO{
   int coef_type; /* all integer, all binary, fractional
			 - considering the type of coefficients*/
   int sign_type; /* same below */
   char var_type; /* '*C'ontinuous,
		     *'B'inary,
		     *'general 'I'nteger,
                     *'F'ixed,
		     *'Z'-continous but can be integerized

		        -those should only appear during preprocessor stage-
                     *negative bina'R'y,
		     *fixable to its 'U'pper bound,
		     *fixable to its 'L'ower bound,
		        -for the last two, need to use is_int to see if
		         they are integer or not-
		     *'T'emporarily fixed,
		     * binary variable and temprarily fixed to
 		       its 'l'ower bound, simiarly,
		       temporarily fixed to its 'u'pper bound
		     */
   int sos_num;      /* #of sos rows that this var appears in */
   int col_size;     /* col size */
   int fix_row_ind; /* state which row caused to fix this variable during
		       basic preprocessor */

   IMPlist *ulist;  /* for binary variables: keeps the list of variables
		       fixed or bounds improved if this variable is fixed to
		       its upper bound */
   IMPlist *llist;  /* same here - lower side */

}COLinfo;

/*===========================================================================*/
/* Data structure to keep relevant info of a row */
/*===========================================================================*/
typedef struct ROWINFO{
   int type; /* all mixed, binary, pure(not binary), cont_binary... */
   int bound_type;  /* all_bounded, mixed
			   - considering the bounds of variables */
   int coef_type; /* all integer, all binary, fractional
			 - considering the type of coefficients*/
   int sign_type; /* all_pos, all_neg, mixed */

   char is_sos_row;
   char * sos_rep;  /* compact representation of the sos row for bitwise
		       operations */

   /* for preprocessor */

   double fixed_obj_offset; /* obtained from fixed vars */
   double fixed_lhs_offset; /* obtained from fixed vars */

   double ub; /* calculated using variable bounds */
   double lb; /* same above */

   double sr_ub; /* calculated using sr relaxations + bounds*/
   double sr_lb; /* same above */

   double orig_ub; /* for debugging purposes */
   double orig_lb;

   int free_var_num;

   int ub_inf_var_num; /* number of variables in this row those cause
			  ub to be infinite */
   int lb_inf_var_num; /* number of variables in this row those cause
			  lb to be infinite */
   int size;
   int fixed_var_num; /* number of fixed variables on this row*/
   int fixable_var_num; /* number of fixable variables on this row*/
   int bin_var_num; /*not fixed binary variables */
   int cont_var_num; /*not fixed continuous variables */
   int frac_coef_num; /* not fixed, frac coeffs on this row */

   char is_redundant;
   char is_updated;
   char vars_checked;

}ROWinfo;

/*===========================================================================*/
/* Data structure to collect information about the model   */
/*===========================================================================*/

typedef struct MIPINFO{
   int prob_type; /* mixed, pure(not binary), binary... */
   int cont_var_num;
   int binary_var_num;
   int binary_var_nz;
   int fixed_var_num;
   int integerizable_var_num;
   int max_row_size;
   int max_col_size;
   int obj_size;  /* number of nonzeros in objective function */

   char is_opt_val_integral; /*is the optimal
			   solution value required to be integral, if one
			   exists*/

   double sum_obj_offset; /* from fixed variables*/

   int binary_sos_row_num; /* sos rows with binary vars count*/
   int binary_row_num; /* rows with binary vars*/
   int cont_row_num; /* rows with cont vars */
   int bin_cont_row_num; /* rows with both cont and bin vars */
   int row_bin_den; /* binary nz / number of rows */
   int col_bin_den; /* binary nz / number of binary columns */
   int row_bin_den_mean; /* 2*row_bin_den*max_row_size/
			    row_bin_den+max_row_size */
   int col_bin_den_mean; /* same here for cols */

   double bin_var_ratio;
   double cont_var_ratio;
   double int_var_ratio;
   double max_row_ratio;
   double max_col_ratio;
   double mat_density;
   double row_density;
   double col_density;
   double sos_bin_row_ratio;
   double bin_row_ratio;

   ROWinfo *rows;
   COLinfo *cols;
}MIPinfo;

/*===========================================================================*/

#if 0
/* not implemented yet */
/* to keep the differences with the original model */
typedef struct MIPDIFF
{
   int rows_del_num;
   int vars_fixed_num;
   int coef_changed_num;
   int bounds_tightened_num;
   int bounds_integerized_num;
   int *rows_deleted_ind;
   int *vars_fixed_ind;
   int *bounds_tightened_ind;
   int *bounds_integerized_ind;
   int *coef_changed_col_ind;
   int *coef_changed_row_ind;
}MIPdiff;

#endif

/*===========================================================================*/
/* This structure stores the user's description of the model */
/*===========================================================================*/
typedef struct MIPDESC{
   int        n;           /* number of columns */
   int        m;           /* number of rows */
   int        nz;          /* number of nonzeros */
   char      *is_int;      /* indicates whether a given variables is integer */
   int       *matbeg;      /* n */
   int       *matind;      /* nz */
   double    *matval;      /* nz */
   double    *obj;         /* n */
   double    *obj1;        /* n */ /* for bicriteria problems */
   double    *obj2;        /* n */ /* for bicriteria problems */
   double    *rhs;         /* m */
   double    *rngval;      /* m */
   char      *sense;       /* m */
   double    *lb;          /* n */
   double    *ub;          /* n */
   char     **colname;     /* column names */
   double     obj_offset;  /* constant to be added to the objective function.*/
   char       obj_sense;   /* objective sense. */

   int        alloc_n;     /* allocated dims */
   int        alloc_m;
   int        alloc_nz;

   int        fixed_n;      /* only used if preprocessor is used */
   int       *fixed_ind;    /* fixed vars to nonzero vals */
   double    *fixed_val;

/* Only to be allocated and used by SYMPHONY */

   int       *col_lengths;
   int       *row_matbeg;      /* m */  /* a row ordered desc for heuristics */
   int       *row_matind;      /* nz */
   double    *row_matval;      /* nz */
   int       *row_lengths;
   /* will keep the orig sense - if prep is used */
   char      *orig_sense;
   int       *orig_ind; /*mapping of indices of presolved model into orig one
			 */

   int        var_type_modified;  /* number of updates on the mip desc */
   int        change_num;  /* number of updates on the mip desc */
   int        change_type[MAX_CHANGE_NUM];  /* type of the mip desc. changes */
   int        new_col_num; /* used only when new cols added */
   int        cru_vars_num;
   int       *cru_vars;
   char       is_modified;

   /* will be evaluated only if preprocessor is used */
   /* it is here to be carried later for further use */
   /* mip info */
   MIPinfo   *mip_inf;

   //  MIPdiff *mip_diff;

}MIPdesc;

/*===========================================================================*\
 * The warm start description contains all information needed to warm start
 * the algorithm.
\*===========================================================================*/

typedef struct WARM_START_DESC{
   bc_node       *rootnode;
   int            cut_num;
   int            allocated_cut_num;
   cut_data     **cuts;
   problem_stat   stat;
   node_times     comp_times;
   int            phase;
   double         lb;
   int            has_ub;
   double         ub;
   lp_sol         best_sol;
   char           trim_tree;
   int            trim_tree_level;
   int            trim_tree_index;
}warm_start_desc;

/*===========================================================================*/
/* solution pool */

typedef struct SP_SOLUTION_DESC{
   double         objval;
   int            xlength;
   int           *xind;
   double        *xval;

   /* The bnb node where this solution was discoverd*/
   int            node_index;

   /* The level of the node in bnb tree where this solution was discovered */
    int            node_level;
}sp_solution;

/*===========================================================================*/

typedef struct SP_DESC{
   /* max. no. of solutions in the pool */
   int            max_solutions;
   /* no. of solutions in the pool */
   int            num_solutions;
   int            total_num_sols_found;
   /* array of those solutions */
   sp_solution    **solutions;
}sp_desc;
#endif
