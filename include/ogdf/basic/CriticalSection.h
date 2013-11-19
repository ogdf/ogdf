/*
 * $Revision: 3840 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-19 08:27:44 +0100 (Di, 19. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of criticial sections.
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

// must be included first here
#include <ogdf/basic/basic.h>


#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_CRITICAL_SECTION_H
#define OGDF_CRITICAL_SECTION_H


#if !defined(OGDF_SYSTEM_WINDOWS)
#include <pthread.h>
#include <errno.h>
#endif

#include <ogdf/basic/System.h>


namespace ogdf {

#if defined(OGDF_SYSTEM_WINDOWS)

//! Representation of a critical section.
/**
 * Critical sections are used to synchronize access to resources shared by
 * several threads. It can be used to protect an object or a piece of code
 * allowing always one thread at a time access.
 */
class OGDF_EXPORT CriticalSection
{
public:
	//! Creates a critical section object.
	CriticalSection() {
#if _WIN32_WINNT >= 0x0600
		InitializeSRWLock(&m_srwLock);
#else
		InitializeCriticalSection(&m_cs);
#endif
	}

	//! Creates a critical section object with spin count.
	/**
	 * The spin count determines how many times the calling thread spins
	 * when the critical section is unavailable, before it performs a wait.
	 * \remark The spin count is only used on multiprocessor systems;
	 *         otherwise it makes no sense.
	 */
	explicit CriticalSection(int spinCount) {
#if _WIN32_WINNT >= 0x0600
#pragma warning( suppress : 4100 )
		InitializeSRWLock(&m_srwLock);
#else
		InitializeCriticalSectionAndSpinCount(&m_cs, spinCount);
#endif
	}

	~CriticalSection() {
#if _WIN32_WINNT < 0x0600
		DeleteCriticalSection(&m_cs);
#endif
	}

	//! Enters the critical section.
	void enter() {
#if _WIN32_WINNT >= 0x0600
		AcquireSRWLockExclusive(&m_srwLock);
#else
		EnterCriticalSection(&m_cs);
#endif
	}

	//! Tries to enter the critical section; returns true on success.
	bool tryEnter() {
#if _WIN32_WINNT >= 0x0600
		return (TryAcquireSRWLockExclusive(&m_srwLock) != 0);
#else
		return (TryEnterCriticalSection(&m_cs) != 0);
#endif
	}

	//! Leaves the critical section.
	void leave() {
#if _WIN32_WINNT >= 0x0600
		ReleaseSRWLockExclusive(&m_srwLock);
#else
		LeaveCriticalSection(&m_cs);
#endif
	}

private:
#if _WIN32_WINNT >= 0x0600
	SRWLOCK m_srwLock;
#else
	CRITICAL_SECTION m_cs; //!< The Windows critical section object.
#endif
};


#else

class OGDF_EXPORT CriticalSection
{
public:
	CriticalSection() : m_spinCount(0) {
		pthread_mutex_init(&m_mutex, NULL);
	}

	explicit CriticalSection(int spinCount) {
		m_spinCount = (System::numberOfProcessors() >= 2) ? spinCount : 0;
		int ret = pthread_mutex_init(&m_mutex, NULL);
		if(ret != 0)
			cout << "initialization of mutex failed: " << ret << endl;
	}

	~CriticalSection() {
		pthread_mutex_destroy(&m_mutex);
	}

	void enter() {
		if(m_spinCount > 0) {
			for(int i = m_spinCount; i > 0; --i)
				if(pthread_mutex_trylock(&m_mutex) != EBUSY)
					return;
		}
		pthread_mutex_lock(&m_mutex);
	}

	bool tryEnter() {
		return (pthread_mutex_trylock(&m_mutex) != EBUSY);
	}

	void leave() {
		pthread_mutex_unlock(&m_mutex);
	}

private:
	pthread_mutex_t m_mutex;
	int             m_spinCount;
};

#endif


} // end namespace ogdf


#endif
