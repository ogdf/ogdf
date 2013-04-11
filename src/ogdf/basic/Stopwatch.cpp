/*
 * $Revision: 3235 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-01-22 15:43:41 +0100 (Di, 22. Jan 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of stopwatch classes
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


#include <ogdf/basic/Stopwatch.h>


namespace ogdf {

	ostream& operator<<(ostream& os, const Stopwatch &stopwatch)
	{
		__int64 centiSeconds = stopwatch.centiSeconds();

		__int64 sec  = centiSeconds/100;
		__int64 mSec = centiSeconds - 100*sec;
		__int64 rSec = sec%60;
		__int64 min  = sec/60;
		__int64 rMin = min%60;

		os << min/60 << ":";
		if(rMin < 10) os << '0';
		os << rMin << ':';
		if(rSec < 10) os << '0';
		os << rSec << '.';
		if (mSec < 10) os << '0';
		os << mSec;
		return os;
	}


	void Stopwatch::start(bool reset)
	{
		if (reset)
			m_totalTime = 0;

		else if (m_running) {
			Logger::ifout() << "Stopwatch::start(): you cannot start a running stopwatch.\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcTimer);
		}

		m_running   = true;
		m_startTime = theTime();
	}


	void Stopwatch::stop()
	{
		if(!m_running) {
			Logger::ifout() << "Stopwatch::stop(): you cannot stop a non-running stopwatch.\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcTimer);
		}

		m_totalTime += theTime() - m_startTime;
		m_running   = false;
	}


	__int64 StopwatchCPU::theTime() const
	{
		double t;
		ogdf::usedTime(t);

		return (__int64)(1000.0*t);
	}


	__int64 StopwatchWallClock::theTime() const
	{
		__int64 t;
		ogdf::System::usedRealTime(t);

		return t;
	}

 }