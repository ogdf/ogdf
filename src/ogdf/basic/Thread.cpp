/*
 * $Revision: 3533 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-06-03 18:22:41 +0200 (Mo, 03. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Thread class representing threads.
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


#include <ogdf/basic/Thread.h>


namespace ogdf {


#ifdef OGDF_SYSTEM_WINDOWS

#ifdef OGDF_USE_THREAD_POOL
	SRWLOCK Thread::s_poolLock = SRWLOCK_INIT;

	int Thread::s_numPoolThreads;
	int Thread::s_numSleepingThreads;
	int Thread::s_maxNumPoolThreads;

	Thread::PoolThreadData **Thread::s_poolThreads;
	Thread::PoolThreadData **Thread::s_sleepingThreads;


	struct Thread::PoolThreadData
	{
		HANDLE              m_handle;		// handle of thread
		unsigned int        m_id;			// thread id
		CONDITION_VARIABLE  m_cvStart;		// signaled when thread shall start executing its work
		Thread             *m_pWork;
	};


	Thread::Thread() : m_handle(0), m_poolThread(0), m_id(0)
	{
		m_evFinished = CreateEvent(NULL, TRUE, FALSE, NULL);
	}


	bool Thread::started() const
	{
		return m_handle || m_poolThread;
	}


	void Thread::start()
	{
		if(started())
			return; // don't start twice

		AcquireSRWLockExclusive(&s_poolLock);

		// sleeping pool thread available?
		if(s_numSleepingThreads > 0) {
			PoolThreadData *pData = s_sleepingThreads[--s_numSleepingThreads];

			m_handle     = 0;
			m_poolThread = pData;
			m_id         = 0;

			pData->m_pWork = this;
			ResetEvent(m_evFinished);
			ReleaseSRWLockExclusive(&s_poolLock);
			WakeConditionVariable(&pData->m_cvStart);

		// create new pool thread?
		} else if(s_numPoolThreads < s_maxNumPoolThreads) {
			PoolThreadData *pData = new PoolThreadData;
			s_poolThreads[s_numPoolThreads++] = pData;

			InitializeConditionVariable(&pData->m_cvStart);

			pData->m_pWork = this;

			m_handle     = 0;
			m_poolThread = pData;
			m_id         = 0;

			pData->m_handle = (HANDLE) _beginthreadex(0, 0, poolThreadProc, pData, 0, &pData->m_id);
			ReleaseSRWLockExclusive(&s_poolLock);

		// create new "normal" thread
		} else {
			ReleaseSRWLockExclusive(&s_poolLock);

			m_poolThread = 0;
			m_handle = (HANDLE) _beginthreadex(0, 0, threadProc, this, 0, &m_id);
		}
	}


	HANDLE Thread::getThreadHandle()
	{
		return m_poolThread ? m_poolThread->m_handle : m_handle;
	}


	long Thread::threadID() const
	{
		return m_poolThread ? (long)m_poolThread->m_id : (long)m_id;
	}


	void Thread::initPool()
	{
		s_maxNumPoolThreads  = System::numberOfProcessors();
		s_numPoolThreads     = 0;
		s_numSleepingThreads = 0;

		s_poolThreads     = new PoolThreadData* [s_maxNumPoolThreads];
		s_sleepingThreads = new PoolThreadData* [s_maxNumPoolThreads];

		for(int i = 0; i < s_maxNumPoolThreads; ++i)
			s_poolThreads[i] = 0;
	}


	void Thread::cleanupPool()
	{
		// test if threads were terminated.
		for (int i = 0; i < s_numPoolThreads; ++i) {
			DWORD exitCode = 0;
			GetExitCodeThread(s_poolThreads[i]->m_handle, &exitCode);
			if (exitCode != STILL_ACTIVE) {
				delete s_poolThreads[i];
				s_poolThreads[i] = 0;
			}
		}

		// wait for all work to be done
		for(int i = 0; i < s_numPoolThreads; ++i) {
			if(s_poolThreads[i] && s_poolThreads[i]->m_pWork != 0) {
				WaitForSingleObject(s_poolThreads[i]->m_pWork->m_evFinished, INFINITE);
			}
		}

		// tell all pool threads to finish
		for(int i = 0; i < s_numPoolThreads; ++i) {
			if (s_poolThreads[i]) {
				WakeConditionVariable(&s_poolThreads[i]->m_cvStart);
			}
		}

		// wait for pool threads to exit
		for(int i = 0; i < s_numPoolThreads; ++i) {
			if (s_poolThreads[i]) {
				WaitForSingleObject(s_poolThreads[i]->m_handle, INFINITE);
			}
		}

		AcquireSRWLockExclusive(&s_poolLock);

		for(int i = 0; i < s_numPoolThreads; ++i) {
			if (s_poolThreads[i]) {
				CloseHandle(s_poolThreads[i]->m_handle);
				delete s_poolThreads[i];
			}
		}

		delete [] s_poolThreads;
		delete [] s_sleepingThreads;

		s_numPoolThreads = 0;

		ReleaseSRWLockExclusive(&s_poolLock);
	}


	unsigned int __stdcall Thread::poolThreadProc(void *pParam)
	{
		PoolThreadData *pData = static_cast<PoolThreadData*>(pParam);

		OGDF_ALLOCATOR::initThread();

		while(pData->m_pWork != 0)
		{
			Thread *pWork = pData->m_pWork;
			pWork->doWork();

			AcquireSRWLockExclusive(&s_poolLock);

			SetEvent(pWork->m_evFinished);
			s_sleepingThreads[s_numSleepingThreads++] = pData;
			pData->m_pWork = 0;

			SleepConditionVariableSRW(&pData->m_cvStart, &s_poolLock, INFINITE, 0);
			ReleaseSRWLockExclusive(&s_poolLock);
		}

		OGDF_ALLOCATOR::flushPool();
		return 0;
	}

#else

	Thread::Thread() : m_handle(0), m_id(0)
	{
		m_evFinished = CreateEvent(NULL, TRUE, FALSE, NULL);
	}


	bool Thread::started() const
	{
		return m_handle != 0;
	}


	void Thread::start()
	{
		if(started())
			return; // don't start twice

		m_handle = (HANDLE) _beginthreadex(0, 0, threadProc, this, 0, &m_id);
	}


	HANDLE Thread::getThreadHandle()
	{
		return m_handle;
	}


	long Thread::threadID() const
	{
		return (long)m_id;
	}


#endif


	Thread::~Thread()
	{
		if(m_handle)
			CloseHandle(m_handle);
		CloseHandle(m_evFinished);
	}


	__uint64 Thread::cpuAffinity(__uint64 mask)
	{
		return SetThreadAffinityMask(getThreadHandle(), (DWORD_PTR)mask);
	}


	void Thread::join()
	{
		WaitForSingleObject(m_evFinished, INFINITE);
	}


	bool Thread::join(unsigned long milliseconds)
	{
		return (WaitForSingleObject(m_evFinished, milliseconds) == WAIT_OBJECT_0);
	}



	unsigned int __stdcall Thread::threadProc(void *pParam)
	{
		Thread *pThread = static_cast<Thread*>(pParam);
		OGDF_ALLOCATOR::initThread();

		pThread->doWork();

		OGDF_ALLOCATOR::flushPool();
		pThread->m_id = 0;

		SetEvent(pThread->m_evFinished);

		return 0;
	}

#else

	Thread::Thread() : m_pt(0) { }

	Thread::~Thread() { }

	bool Thread::started() const {
		return m_pt != 0;
	}


	long Thread::threadID() const
	{
		return (long)m_pt;
	}


	// not supported
	__uint64 Thread::cpuAffinity(__uint64 mask) { return mask; }

	void Thread::start()
	{
		OGDF_ASSERT(m_pt == 0);
		pthread_create(&m_pt, NULL, threadProc, this);
	}


	void Thread::join()
	{
		if(m_pt != 0)
			pthread_join(m_pt,NULL);
	}


	void *Thread::threadProc(void *pParam)
	{
		Thread *pThread = static_cast<Thread*>(pParam);
		OGDF_ALLOCATOR::initThread();
		pThread->doWork();
		//pthread_exit(NULL);
		OGDF_ALLOCATOR::flushPool();
		pThread->m_pt = 0;
		return 0;
	}

#endif


}
