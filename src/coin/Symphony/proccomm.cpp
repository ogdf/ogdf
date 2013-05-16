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

#include <stdio.h>
#include <stdlib.h>

#include "sym_proccomm.h"

#ifdef __PVM__

#include <pvm3.h>

/*===========================================================================*
 * These functions are just a front-end to PVM
 *===========================================================================*/

/*===========================================================================*/

int register_process()
{
   int mytid;

   PVM_FUNC(mytid, pvm_mytid());

   return(mytid);
}

/*===========================================================================*/

int init_send(int data_packing)
{
   int s_bufid;

   PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));

   return(s_bufid);
}

/*===========================================================================*/

int send_char_array(char *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_pkbyte(array, size, 1));

   return(info);
}

/*===========================================================================*/

int send_int_array(int *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_pkint(array, size, 1));

   return(info);
}

/*===========================================================================*/

int send_float_array(float *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_pkfloat(array, size, 1));

   return(info);
}

/*===========================================================================*/

int send_dbl_array(double *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_pkdouble(array, size, 1));

   return(info);
}

/*===========================================================================*/

int send_str(char *str)
{
   int info;

   PVM_FUNC(info, pvm_pkstr(str));

   return(info);
}

/*===========================================================================*/

int send_msg(int recipient, int msgtag)
{
   int info;

   PVM_FUNC(info, pvm_send(recipient, msgtag));

   return(info);
}

/*===========================================================================*/

int msend_msg(int *recipients, int number, int msgtag)
{
   int info;

   PVM_FUNC(info, pvm_mcast(recipients, number, msgtag));

   return(info);
}

/*===========================================================================*/

int receive_msg(int who, int what)
{
   int r_bufid;

   PVM_FUNC(r_bufid, pvm_recv(who, what));

   return(r_bufid);
}

/*===========================================================================*/

int nreceive_msg(int who, int what)
{
   int r_bufid;

   PVM_FUNC(r_bufid, pvm_nrecv(who, what));

   return(r_bufid);
}

/*===========================================================================*/

int treceive_msg(int who, int what, struct timeval *timeout)
{
   int r_bufid;

   PVM_FUNC(r_bufid, pvm_trecv(who, what, timeout));

   return(r_bufid);
}

/*===========================================================================*/

int bufinfo(int r_bufid, int *bytes, int *msgtag, int *sender)
{
   int info;

   PVM_FUNC(info, pvm_bufinfo(r_bufid, bytes, msgtag, sender));

   return(info);
}

/*===========================================================================*/

int freebuf(int bufid)
{
   int info;

   PVM_FUNC(info, pvm_freebuf(bufid));

   return(info);
}

/*===========================================================================*/

int receive_char_array(char *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_upkbyte(array, size, 1));

   return(info);
}

/*===========================================================================*/

int receive_int_array(int *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_upkint(array, size, 1));

   return(info);
}

/*===========================================================================*/

int receive_dbl_array(double *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_upkdouble(array, size, 1));

   return(info);
}

/*===========================================================================*/

int receive_float_array(float *array, int size)
{
   int info;

   PVM_FUNC(info, pvm_upkfloat(array, size, 1));

   return(info);
}

/*===========================================================================*/

int receive_str(char *str)
{
   int info;

   PVM_FUNC(info, pvm_upkstr(str));

   return(info);
}

/*===========================================================================*/

void comm_exit(void)
{
   pvm_exit();
}

/*===========================================================================*/

int spawn(char *task, char **argv, int flag, char *where, int ntask,
	  int *tids)
{
   int status;

   if ((status = pvm_spawn(task, argv, flag, where, ntask, tids)) != ntask){
      printf("Couldn't start %s!! \n", task);
      PVM_ERROR(*tids);
      exit(0);
   }
   return(status);
}

/*===========================================================================*/

int pstat(int tid)
{
   return(pvm_pstat(tid));
}

/*===========================================================================*/

void kill_proc(int tid)
{
   pvm_kill(tid);
}

/*===========================================================================*/

void setsbuf(int sbufid)
{
   int info;

   PVM_FUNC(info, pvm_setsbuf(sbufid));
}

/*===========================================================================*/

void setrbuf(int rbufid)
{
   int info;

   PVM_FUNC(info, pvm_setrbuf(rbufid));
}

/*===========================================================================*/

void PVM_ERROR(int info)
{
   printf("Pvm Error %i : ", info);
   switch (info){
      case PvmOk: printf("PvmOk(Error 0)"); break;
      case PvmBadParam: printf("PvmBadParam (Bad parameter)"); break;
      case PvmMismatch: printf("PvmMismatch (Count mismatch)"); break;
      case PvmOverflow: printf("PvmOverflow (Value too large)"); break;
      case PvmNoData: printf("PvmNoData (End of buffer)"); break;
      case PvmNoHost: printf("PvmNoHost (No such host)"); break;
      case PvmNoFile: printf("PvmNoFile (No such file)"); break;
      case PvmNoMem: printf("PvmNoMem (Malloc failed)"); break;
      case PvmBadMsg: printf("PvmBadMsg (Can't decode message)"); break;
      case PvmSysErr: printf("PvmSysErr (Can't contact local daemon)"); break;
      case PvmNoBuf: printf("PvmNoBuf (No current buffer)"); break;
      case PvmNoSuchBuf: printf("PvmNoSuchBuf (No such buffer)"); break;
      case PvmNullGroup: printf("PvmNullGroup (Null group name)"); break;
      case PvmDupGroup: printf("PvmDupGroup (Already in group)"); break;
      case PvmNoGroup: printf("PvmNoGroup (No such group)"); break;
      case PvmNotInGroup: printf("PvmNotInGroup (Not in group)"); break;
      case PvmNoInst: printf("PvmNoInst (No such instance)"); break;
      case PvmHostFail: printf("PvmHostFail (Host failed) "); break;
      case PvmNoParent: printf("PvmNoParent (No parent task)"); break;
      case PvmNotImpl: printf("PvmNotImpl (Not implemented)"); break;
      case PvmDSysErr: printf("PvmDSysErr (Pvmd system error)"); break;
      case PvmBadVersion: printf("PvmBadVersion (Version mismatch)"); break;
      case PvmOutOfRes: printf("PvmOutOfRes (Out of resources)"); break;
      case PvmDupHost: printf("PvmDupHost (Duplicate host)"); break;
      case PvmCantStart: printf("PvmCantStart (Can't start pvmd)"); break;
      case PvmAlready: printf("PvmAlready (Already in progress)"); break;
      case PvmNoTask: printf("PvmNoTask (No such task)"); break;
#if 0
      case PvmNoEntry: printf("PvmNoEntry (No such entry)"); break;
      case PvmDupEntry: printf("PvmDupEntry (Duplicate entry)"); break;
#endif
   }
   printf("\n\n");
   fflush(stdout);
}

#else

/*===========================================================================*
 * Function stubs in case PVM is not used
 *===========================================================================*/

/*===========================================================================*/

int register_process()
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int init_send(int data_packing)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_char_array(char *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_int_array(int *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_dbl_array(double *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_float_array(float *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_str(char *str)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int send_msg(int sender, int msgtag)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int msend_msg(int *recipients, int number, int msgtag)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}
/*===========================================================================*/

int receive_msg(int who, int what)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int nreceive_msg(int who, int what)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int treceive_msg(int who, int what, struct timeval *timeout)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int bufinfo(int r_bufid, int *bytes, int *msgtag, int *sender)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int freebuf(int bufid)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int receive_char_array(char *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int receive_int_array(int *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int receive_dbl_array(double *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int receive_float_array(float *array, int size)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int receive_str(char *str)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

void comm_exit(void)
{
   printf("\nComm Error: Unknown communications protocol\n\n");
}

/*===========================================================================*/

int spawn(char *task, char **argv, int flag, char *where, int ntask,
	  int *tids)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

int pstat(int tid)
{
   printf("\nComm Error: Unknown communications protocol\n\n");

   return(0);
}

/*===========================================================================*/

void kill_proc(int tid)
{
   printf("\nComm Error: Unknown communications protocol\n\n");
}

/*===========================================================================*/

void setsbuf(int sbufid)
{
   printf("\nComm Error: Unknown communications protocol\n\n");
}

/*===========================================================================*/

void setrbuf(int rbufid)
{
   printf("\nComm Error: Unknown communications protocol\n\n");
}

/*===========================================================================*/

#endif
