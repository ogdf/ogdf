/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of Thread class representing threads.
 *
 * \author Carsten Gutwenger
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_THREAD_H
#define OGDF_THREAD_H

#include <ogdf/basic/basic.h>


#ifdef OGDF_SYSTEM_WINDOWS
#include <process.h>

#if _WIN32_WINNT >= 0x0600
#define OGDF_USE_THREAD_POOL
#endif

#else
#include <pthread.h>
#endif


namespace ogdf {


//! Base class for threads.
class Thread
{
	friend class Initialization;

#ifdef OGDF_USE_THREAD_POOL
	static void initPool();
	static void cleanupPool();
#endif

public:

	//! Initializes a thread, but does not create a system thread yet.
	Thread();

	//! Destructor. Frees resources.
	virtual ~Thread();

	//! Returns the ID of the system thread associated with the thread object.
	long threadID() const;

	//! Returns whether the thread has been started (i.e. a system thread is currently associated with the thread object).
	bool started() const;

	//! Sets the CPU affinity mask of the thread to \a mask.
	__uint64 cpuAffinity(__uint64 mask);

	//! Starts execution of the thread.
	void start();

	//! Waits until the thread has finished.
	void join();

	//! Waits until the thread has finished or the time-out interval of \a milliseconds elapses.
	/**
	 * @param milliseconds is the time-out interval in milliseconds.
	 * @return true if the thread has finished or false if the time-out interval has elapsed.
	 */
	bool join(unsigned long milliseconds);

protected:
	//! The actual work perfomed by the thread. Must be defined by derived classes.
	virtual void doWork() = 0;

private:

#ifdef OGDF_SYSTEM_WINDOWS

#ifdef OGDF_USE_THREAD_POOL
	struct PoolThreadData;

	static SRWLOCK s_poolLock;
	static int s_numPoolThreads;
	static int s_numSleepingThreads;
	static int s_maxNumPoolThreads;

	static PoolThreadData **s_poolThreads;
	static PoolThreadData **s_sleepingThreads;

	PoolThreadData *m_poolThread;	//!< associated pool thread (0 if no pool thread is used)
#endif

	HANDLE          m_handle;		//!< thread handle (0 if pool thread)
	unsigned int    m_id;			//!< thread id (0 if pool thread)
	HANDLE          m_evFinished;	//!< event which is signaled by thread when work is finished

	HANDLE getThreadHandle();

	//! Thread procedure for normal threads; \a pParam points to the Thread object.
	static unsigned int __stdcall threadProc(void *pParam);

#ifdef OGDF_USE_THREAD_POOL
	//! Thread procedure for pool threads; \a pParam points to the PoolThreadData object.
	static unsigned int __stdcall poolThreadProc(void *pParam);
#endif


#else

	static void *threadProc(void *pParam);

	pthread_t m_pt;

#endif

	OGDF_NEW_DELETE

};


} // end namespace ogdf


#endif
