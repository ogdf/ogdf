/** \file
 * \brief Simple, safe base classes for C++ observables and observers.
 *
 * \author Simon D. Fink
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#include <ogdf/basic/List.h>
#include <ogdf/basic/internal/config.h>

#ifndef OGDF_MEMORY_POOL_NTS
#	include <mutex>
#endif

namespace ogdf {

/**
 * Base class for an observer for a single Observable object.
 * Will be notified when the observable is destructed and can be subclassed to provide further callbacks.
 * For compatibility with MSVC, the Observer subclass has to be defined before the Observable subclass.
 *
 * @tparam TObserved The subclass of Observable that will be observed.
 * @tparam TObserver The subclass of Observer that defines all virtual callback methods.
 */
template<typename TObserved, typename TObserver>
class Observer {
public:
	//! Constructs instance of Observer class
	Observer() { }

	//! Destroys the instance, unregisters it from watched instance.
	virtual ~Observer() { reregister(nullptr); }

	//! Associates observer instance with instance \p obs.
	void reregister(const TObserved* obs) {
		const TObserved* old = m_pObserved;
		if (m_pObserved) {
			m_pObserved->unregisterObserver(m_itObsList);
		}
		m_pObserved = obs;
		if (m_pObserved != nullptr) {
			m_itObsList = obs->registerObserver(dynamic_cast<TObserver*>(this));
		}
		registrationChanged(old);
	}

	//! Called after reregister() changed the observed instance.
	virtual void registrationChanged(const TObserved* old) { }

	const TObserved* getObserved() const { return m_pObserved; }

private:
	const TObserved* m_pObserved = nullptr; //! watched instance
	typename ListPure<TObserver*>::iterator m_itObsList; //! own entry in m_pObserved's observer list
};

/**
 * Base class for an observable object that can be tracked by multiple Observer objects.
 * Will be notify its observers when it is destructed and can be subclassed to provide further callbacks.
 * For compatibility with MSVC, the Observer subclass has to be defined before the Observable subclass.
 *
 * @tparam TObserved The subclass of Observable that will be observed.
 * @tparam TObserver The subclass of Observer that defines all virtual callback methods.
 */
template<typename TObserver, typename TObserved>
class Observable {
	friend Observer<TObserved, TObserver>;

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared acces to register/unregister methods.
#endif
	mutable ListPure<TObserver*> m_regObservers; //!< The registered observers.

public:
	virtual ~Observable() {
		// clearObservers must be called by child class, calling it here would
		// notify Observers with an already partially destructed Observable
		OGDF_ASSERT(m_regObservers.empty());
	}

private:
	//! Registers an observer.
	/**
	 * @param obs is a pointer to the observer that shall be registered
	 * @return an iterator pointing to the entry for the registered observer in the list of registered
	 *         observers. This iterator is required for unregistering the observer again.
	 */
	typename ListPure<TObserver*>::iterator registerObserver(TObserver* obs) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		OGDF_ASSERT(obs != nullptr);
		return m_regObservers.pushBack(obs);
	}

	//! Unregisters an observer.
	/**
	 * @param it is an iterator pointing to the entry in the list of registered observers for the
	 *           observer to be unregistered.
	 */
	void unregisterObserver(typename ListPure<TObserver*>::iterator it) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_regObservers.del(it);
	}

protected:
	const ListPure<TObserver*>& getObservers() const { return m_regObservers; }

	void clearObservers() {
		while (!m_regObservers.empty()) {
			TObserver* obs = m_regObservers.front();
			obs->reregister(nullptr);
			OGDF_ASSERT(m_regObservers.empty()
					|| (m_regObservers.front() != obs && m_regObservers.back() != obs));
		}
	}
};
}
