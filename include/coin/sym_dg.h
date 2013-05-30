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

#ifndef _INTERMEDIARY_H
#define _INTERMEDIARY_H

#include <stdio.h>

#include "sym_proto.h"
#include "sym_dg_params.h"
#include "sym_dg_u.h"

#ifndef PIPE_ERROR
#define PIPE_ERROR {printf("PIPE_ERROR !!!!!!!!!!!"); }
#endif

/*===========================================================================*/

typedef struct WINDOW_DESCRIPTION{
   int         canvas_width;
   int         canvas_height;
   int         viewable_width;
   int         viewable_height;
   int         disp_nodelabels;
   int         disp_nodeweights;
   int         disp_edgeweights;
   char        node_dash[MAX_DASH_PATTERN_LENGTH +1];
   char        edge_dash[MAX_DASH_PATTERN_LENGTH +1];
   int         node_radius;
   int         interactive_mode;
   int         mouse_tracking;
   double      scale_factor;
   char        nodelabel_font[MAX_FONT_LENGTH +1];
   char        nodeweight_font[MAX_FONT_LENGTH +1];
   char        edgeweight_font[MAX_FONT_LENGTH +1];
}win_desc;

/*===========================================================================*/

typedef struct DG_NODE{
   int     node_id;
   int     posx;
   int     posy;
   int     radius;
   char    deleted;   /* 1/0 */
   char    label[MAX_LABEL_LENGTH +1]; /*watch out, terminating character
				      is a zero!*/
   char    weight[MAX_WEIGHT_LENGTH +1];
   char    dash[MAX_DASH_PATTERN_LENGTH +1];
}dg_node;

/*===========================================================================*/

typedef struct DG_EDGE{
   int     edge_id;
   int     tail;
   int     head;
   char    deleted;
   char    weight[MAX_WEIGHT_LENGTH +1];
   char    dash[MAX_DASH_PATTERN_LENGTH +1];
}dg_edge;

/*===========================================================================*/

typedef struct DG_GRAPH{
   int         nodenum;
   int         deleted_nodenum;
   dg_node    *nodes;
   int         edgenum;
   int         deleted_edgenum;
   dg_edge    *edges;
}dg_graph;

/*===========================================================================*/

typedef struct BUF_FIFO{
   int        *bufid;
   int         bufspace;
   int         bufwrite;
   int         bufread;
}buf_fifo;

/*===========================================================================*/

typedef struct WINDOW{
   char        name[MAX_NAME_LENGTH +1];
   unsigned int id;
   void       *user;

   int         owner_tid;     /* tid of process that initiated this */
   buf_fifo    buf;      /* structure to store buffer id's in a FIFO */
   int         window_displayed;  /* T/F T if window is displayed */
   int         wait_for_click;  /* 0:don't wait, 1:wait, no report,
				   2:wait and report */
   win_desc    desc;
   char        title[MAX_TITLE_LENGTH +1];
   dg_graph    g;
   int         copy_status; /*0:nothing, 1:waiting to be copied, 2:waiting to
			      get a copy. */
   char        source[MAX_NAME_LENGTH +1]; /* set if copy_status=2 */
   char        target[MAX_NAME_LENGTH +1]; /* set if copy_status=1 */

   int         text_length;
   char       *text;
}window;

/*===========================================================================*/

typedef struct DRAW_GRAPH_PROBLEM{
   void         *user;
   int           master;
   dg_params     par;
   char          waiting_to_die;
   unsigned int  next_id;
   int           window_num;    /* the number of windows */
   window      **windows;

   int           message_length;
   char         *message;
}dg_prob;

/*===========================================================================*/

void INTERMED_ERROR PROTO((char *window_name, int old_msgtag,
			   int receiver, int msgtag));
int spprint PROTO((FILE *write_to, const char *format, ...));

int start_child PROTO((char *cmd, FILE **readpipe, FILE **writepipe));
int find_window PROTO((int window_num, window **windows, char *name));
void read_node_desc_from_pvm PROTO((dg_node *nod, window *win));
void read_edge_desc_from_pvm PROTO((dg_edge *edg, window *win));
int find_node PROTO((int node_id, dg_graph *g));
int find_edge PROTO((int edge_id, dg_graph *g));
void compress_graph PROTO((dg_graph *g));
void copy_window_structure PROTO((window *target_win, window *source_win));
void display_graph_on_canvas PROTO((window *win, FILE *write_to));
void free_window PROTO((int *pwindow_num, window **windows, int i));
void copy_win_desc_from_par PROTO((window *win, dg_params *par));
void set_window_desc_pvm PROTO((int key, window *win));
void wait_for_you_can_die PROTO((dg_prob *dgp, FILE *write_to));
window *init_dgwin PROTO((dg_prob *dgp, int sender, char *name, char *title));
void add_msg PROTO((window *win, int bufid));
int get_next_msg PROTO((window *win));

//FILE *fdopen PROTO((int, const char *));

#endif
