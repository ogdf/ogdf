/*
 * $Revision: 3257 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-01-25 10:41:13 +0100 (Fr, 25. Jan 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of stopwatch classes
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

#ifndef OGDF_STOPWATCH_H
#define OGDF_STOPWATCH_H


#include <ogdf/basic/basic.h>


namespace ogdf {

	class OGDF_EXPORT Stopwatch {

		__int64 m_startTime;  //!< The start time of the timer in milliseconds.
		__int64 m_totalTime;  //!< The total time in milliseconds.
		bool    m_running;    //!< true, if the timer is running.

	public:
		//! Initializes a stop watch with total time 0.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 */
		Stopwatch() : m_totalTime(0), m_running(false) { }

		//! Initializes a stopwatch and sets its total time to \a milliSecs.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 *
		 * \param milliSecs The intial value of the total time in milliseconds.
		 */
		Stopwatch(__int64 milliSecs) : m_totalTime(milliSecs), m_running(false) { }


		virtual ~Stopwatch() { }


		//! Starts the stopwatch.
		/**
		 * For safety reasons starting a running timer is an error.
		 *
		 * \param reset If this flag is set to true, the stopwatch is reset before it is started.
		 */
		void start(bool reset = false);

		//! Stops the stopwatch and adds the difference between the current time and the starting time to the total time.
		/**
		 * Stopping a non-running stopwatch is an error.
		 */
		void stop();

		//! Stops the stopwatch and sets its total time to 0.
		void reset() {
			m_running   = false;
			m_totalTime = 0;
		}


		//! Returns true if the stopwatch is running, false otherwise.
		bool running() const { return m_running; }


		//! Returns the currently elapsed time in milliseconds.
		/**
		 * It is not necessary to stop the timer to get the correct time.
		 */
		__int64 milliSeconds() const {
			return (m_running) ? (m_totalTime + theTime() - m_startTime) : m_totalTime;
		}

		//! Returns the currently elapsed time in 1/100-seconds.
		/**
		 * It is not necessary to stop the timer to get the correct time.
		 */
		__int64 centiSeconds() const { return milliSeconds()/10; }

		//! Returns the currently elapsed time in seconds.
		/**
		 * It is not necessary to stop the timer to get the correct time.
		 * The result is rounded down to the next integer value.
		 */
		__int64 seconds() const { return milliSeconds()/1000; }

		//! Returns the currently elapsed time in minutes.
		/**
		 * It is not necessary to stop the timer to get the correct time.
		 * The result is rounded down to the next integer value.
		 */
		__int64 minutes() const { return seconds()/60; }

		//! Returns the currently elapsed time in hours.
		/**
		 * It is not necessary to stop the timer to get the correct time.
		 * The result is rounded down to the next integer value.
		 */
		__int64 hours() const { return seconds()/3600; }


		//! Returns true iff the currently elapsed time exceeds \a maxSeconds.
		bool exceeds(__int64 maxSeconds) const {
			return (seconds() >= maxSeconds);
		}

		//! Adds \a centiSeconds to total time.
		/**
		 * \param centiSeconds The number of centiseconds to be added.
		 */
		void addCentiSeconds(__int64 centiSeconds) {
			m_totalTime += 10*centiSeconds;
		}


		//! Writes the currently elapsed time in the format <tt>hh:mm:ss.sec/100</tt> to output stream \a os.
		/**
		 * \param os        The output stream.
		 * \param stopwatch The stopwatch whose elapsed time shall be written.
		 * \return A reference to the output stream \a os.
		 */
		friend ostream& operator<<(ostream& os, const Stopwatch &stopwatch);

	protected:

		//! Returns the current time in milliseconds (from some fixed starting point).
		/**
		 * This pure virtual function is used for measuring time differences.
		 * It must be implemented in derived classes.
		 *
		 * \return The time since some starting point (e.g., the program start) in milliseconds.
		 */
		virtual __int64 theTime() const = 0;

	};


	//! Implements a stopwatch measuring CPU time.
	class OGDF_EXPORT StopwatchCPU : public Stopwatch {
	public:
		//! Creates a stopwatch for measuring CPU time with total time 0.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 */
		StopwatchCPU() : Stopwatch() { }

		//! Creates a stopwatch for measuring CPU time and sets its total time to \a milliSecs.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 *
		 * \param milliSecs The intial value of the total time in milliseconds.
		 */
		StopwatchCPU(__int64 milliSecs) : Stopwatch(milliSecs) { }

		virtual ~StopwatchCPU() { }

	private:
		//! Returns the current CPU time in milliseconds (from some fixed starting point).
		virtual __int64 theTime() const;
	};


	//! Implements a stopwatch measuring wall-clock time.
	class OGDF_EXPORT StopwatchWallClock : public Stopwatch {
	public:
		//! Creates a stopwatch for measuring wall-clock time with total time 0.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 */
		StopwatchWallClock() : Stopwatch() { }

		//! Creates a stopwatch for measuring wall-clock time and sets its total time to \a milliSecs.
		/**
		 * After creation the stopwatch is not running, i.e., it has to be started explicitly
		 * for measuring time.
		 *
		 * \param milliSecs The intial value of the total time in milliseconds.
		 */
		StopwatchWallClock(__int64 milliSecs) : Stopwatch(milliSecs) { }

		virtual ~StopwatchWallClock() { }

	private:
		//! Returns the current wall-clock time in milliseconds (from some fixed starting point).
		virtual __int64 theTime() const;
	};

}


#endif
