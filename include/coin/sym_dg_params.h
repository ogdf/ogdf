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

#ifndef _DRAW_GRAPH_PARAMS_H
#define _DRAW_GRAPH_PARAMS_H

#include "sym_proto.h"

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************         Constants related to DrawGraph               **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

#define MAX_HORIZ                         500
#define MAX_VERT                          500
#define RADIUS                            10

#define  DONT_WAIT_FOR_CLICK              0
#define  WAIT_FOR_CLICK_AND_REPORT        1
#define  WAIT_FOR_CLICK_NO_REPORT         2


#define MAX_NAME_LENGTH                   36
#define MAX_TITLE_LENGTH                  64
#define MAX_DASH_PATTERN_LENGTH           32
#define MAX_FONT_LENGTH                   128
#define MAX_LABEL_LENGTH                  8
#define MAX_WEIGHT_LENGTH                 12


#define IGDTOI_CLICK_HAPPENED             90100
#define IGDTOI_QUIT_WINDOW                90101
#define IGDTOI_QUIT_APPLICATION           90102
#define IGDTOI_TEXT_ENTERED               90103
#define IGDTOI_REQUEST_GRAPH              90104


#define ITOC_CLICK_HAPPENED               90200
#define ITOC_WINDOW_DOESNT_EXIST          90201
#define ITOC_WINDOW_ALREADY_EXISTS        90202
#define ITOC_WINDOW_ISNT_DISPLAYED        90203
#define ITOC_APPLICATION_KILLED           90204
#define ITOC_USER_MESSAGE                 90205


#define CTOI_INITIALIZE_WINDOW            90300
#define CTOI_QUIT_WINDOW                  90301
#define CTOI_QUIT                         90302
#define CTOI_YOU_CAN_DIE                  90303
#define CTOI_COPY_GRAPH                   90304
#define CTOI_CHANGE_WINDOW_DESC           90305
#define CTOI_SET_GRAPH                    90306
#define CTOI_DRAW_GRAPH                   90307
#define CTOI_SET_AND_DRAW_GRAPH           90308
#define CTOI_DELETE_GRAPH                 90309
#define CTOI_WAIT_FOR_CLICK_NO_REPORT     90310
#define CTOI_WAIT_FOR_CLICK_AND_REPORT    90311
#define CTOI_SAVE_GRAPH_TO_FILE           90312
#define CTOI_SAVE_GRAPH_PS_TO_FILE        90313
#define CTOI_CLONE_WINDOW                 90314
#define CTOI_RENAME_WINDOW                90315
#define CTOI_RESIZE_VIEWABLE_WINDOW       90316
#define CTOI_RESIZE_CANVAS                90317
#define CTOI_MODIFY_GRAPH                 90318
#define CTOI_CLEAR_MESSAGE                90319
#define CTOI_PRINT_MESSAGE                90320
#define CTOI_APPEND_MESSAGE               90321

#define CTOI_USER_MESSAGE                 90350


#define WAITING_TO_GET_A_COPY             90400
#define WAITING_TO_BE_COPIED              90401
#define MODIFY_ADD_NODES                  90402
#define MODIFY_CHANGE_WEIGHTS_OF_NODES    90403
#define MODIFY_CHANGE_LABELS_OF_NODES     90404
#define MODIFY_CHANGE_DASH_OF_NODES       90405
#define MODIFY_CHANGE_RADII_OF_NODES      90406
#define MODIFY_DELETE_NODES               90407
#define MODIFY_ADD_EDGES                  90408
#define MODIFY_CHANGE_WEIGHTS_OF_EDGES    90409
#define MODIFY_CHANGE_DASH_OF_EDGES       90410
#define MODIFY_DELETE_EDGES               90411
#define MODIFY_DELETE_ALL_EDGES           90412
#define MODIFY_END_OF_MESSAGE             90413

/* Window description keycodes */
#define CANVAS_WIDTH                      1
#define CANVAS_HEIGHT                     2
#define VIEWABLE_WIDTH                    3
#define VIEWABLE_HEIGHT                   4
#define DISP_NODELABELS                   5
#define DISP_NODEWEIGHTS                  6
#define DISP_EDGEWEIGHTS                  7
#define NODE_DASH                         8
#define EDGE_DASH                         9
#define NODE_RADIUS                       10
#define INTERACTIVE_MODE                  11
#define MOUSE_TRACKING                    12
#define SCALE_FACTOR                      13
#define NODELABEL_FONT                    14
#define NODEWEIGHT_FONT                   15
#define EDGEWEIGHT_FONT                   16

typedef struct DRAW_GRAPH_PARAMS{
   char   source_path[MAX_FILE_NAME_LENGTH +1];
   int    echo_commands;

   int    canvas_width;
   int    canvas_height;
   int    viewable_width;
   int    viewable_height;
   int    disp_nodelabels;
   int    disp_nodeweights;
   int    disp_edgeweights;
   char   node_dash[MAX_DASH_PATTERN_LENGTH +1];
   char   edge_dash[MAX_DASH_PATTERN_LENGTH +1];
   int    node_radius;
   int    interactive_mode;
   int    mouse_tracking;
   double scale_factor;
   char   nodelabel_font[MAX_FONT_LENGTH +1];
   char   nodeweight_font[MAX_FONT_LENGTH +1];
   char   edgeweight_font[MAX_FONT_LENGTH +1];
}dg_params;

#endif
