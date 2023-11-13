#pragma once

#include <ogdf/basic/List.h>
#include <ogdf/basic/internal/config.h>

#ifndef OGDF_MEMORY_POOL_NTS
#	include <mutex>
#endif

namespace ogdf {

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

	//! Called after ::reregister changed the observed instance.
	virtual void registrationChanged(const TObserved* old) { }

	const TObserved* getObserved() const { return m_pObserved; }

private:
	const TObserved* m_pObserved = nullptr; //! watched instance
	typename ListPure<TObserver*>::iterator m_itObsList; //! own entry in m_pObserved's observer list
};

template<typename TObserver, typename TObserved>
class Observable {
	friend Observer<TObserved, TObserver>;

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared acces to register/unregister methods.
#endif
	mutable ListPure<TObserver*> m_regObservers; //!< The registered observers.

public:
	virtual ~Observable() {
		// clearAllObservers must be called by child class, calling it here would
		// notify Observers with an already partially destructed Observable
		OGDF_ASSERT(m_regObservers.empty());
	}

private:
	//! Registers a observer.
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

	//! Unregisters a observer.
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
