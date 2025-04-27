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
#include <ogdf/basic/basic.h>
#include <ogdf/basic/internal/copy_move.h>

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
	//! Constructs unregistered instance of Observer class
	/**
	 * Note that calling reregister in the constructor of some intermediate class will trigger
	 * registrationChanged while child classes are not yet fully constructed.
	 */
	Observer() = default;

	OGDF_DEPRECATED("calls registrationChanged with only partially-constructed child classes, "
					"see copy constructor of Observer for fix")

	explicit Observer(const TObserved* R) { reregister(R); }

	//! Destroys the instance, unregisters it from watched instance.
	/**
	 * Callback registrationChanged will not be made from destructor as all subclasses are already
	 * partially destroyed at that point.
	 */
	virtual ~Observer() {
		if (m_pObserved) {
			m_pObserved->unregisterObserver(m_itObsList);
			m_pObserved = nullptr;
		}
	}

	/**
	 * If you want to copy a subclass of Observer, call the default Observer() constructor and
	 * optionally also call reregister if it makes sense.
	 */
	OGDF_NO_COPY(Observer)
	/**
	 * If you want to move a subclass of Observer, call the default Observer() constructor and
	 * optionally also call reregister if it makes sense.
	 */
	OGDF_NO_MOVE(Observer)

	//! Associates observer instance with instance \p obs.
	/**
	 * This always unregisters and reregisters the observer, even if \p obs == getObserved().
	 * Consequently, this observer will always be the last in the list to be notified of events.
	 * Furthermore, registrationChanged() will always be fired when this method is called.
	 */
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
 * Will notify its observers when it is destructed and can be subclassed to provide further callbacks.
 * For compatibility with MSVC, the Observer subclass has to be defined before the Observable subclass.
 *
 * @tparam TObserved The subclass of Observable that will be observed.
 * @tparam TObserver The subclass of Observer that defines all virtual callback methods.
 */
template<typename TObserver, typename TObserved>
class Observable {
	friend Observer<TObserved, TObserver>;
	friend TObserved;
	friend TObserver;

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared access to register/unregister methods.
#endif
	mutable ListPure<TObserver*> m_regObservers; //!< The registered observers.

public:
	Observable() = default;

	/**
	 * Note that all Observers must already be removed once the destructor of this base class is
	 * invoked (e.g. through clearObservers in a child class destructor), as calling it here would
	 * notify Observers with already partially destructed child classes.
	 */
	virtual ~Observable() { OGDF_ASSERT(m_regObservers.empty()); }

	/**
	 * If you want to copy a subclass of Observable, call the default Observable() constructor.
	 * Note that Observers can only observe one Observable and in this case will stay with the old one.
	 */
	OGDF_NO_COPY(Observable)
	/**
	 * If you want to move a subclass of Observable, call the default Observable() constructor.
	 * If you want to also point all Observers to the new location of their Observable, call
	 * register for each of them.
	 */
	OGDF_NO_MOVE(Observable)

protected:
	//! Registers an observer.
	/**
	 * You should never directly call this method as it is called by the Observer() constructor and
	 * Observer::reregister() methods automatically.
	 * If you have a class that inherits from multiple Observables, you may need to reexport this method:
	 * \code
	 * friend Observer<ObservableA, ObserverA>;
	 * friend Observer<ObservableB, ObserverB>;
	 * using ParentObservableA::(un)registerObserver; using ParentObservableB::(un)registerObserver;
	 * \endcode
	 * See https://stackoverflow.com/a/1313162 and ClusterGraph for an example.
	 *
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
	 * You should never directly call this method as it is called by the ~Observer() destructor and
	 * Observer::reregister() methods automatically.
	 * If you have a class that inherits from multiple Observables, you may need to reexport this method:
	 * \code
	 * friend Observer<ObservableA, ObserverA>;
	 * friend Observer<ObservableB, ObserverB>;
	 * using ParentObservableA::(un)registerObserver; using ParentObservableB::(un)registerObserver;
	 * \endcode
	 * See https://stackoverflow.com/a/1313162 and ClusterGraph for an example.
	 *
	 * @param it is an iterator pointing to the entry in the list of registered observers for the
	 *           observer to be unregistered.
	 */
	void unregisterObserver(typename ListPure<TObserver*>::iterator it) const {
#ifndef OGDF_MEMORY_POOL_NTS
		std::lock_guard<std::mutex> guard(m_mutexRegArrays);
#endif
		m_regObservers.del(it);
	}

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
