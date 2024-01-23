#pragma once

#include <ogdf/basic/pctree/PCEnum.h>

namespace pc_tree {
class PCTreeForest;

template<class Key>
class PCTreeRegistry : public ogdf::RegistryBase<Key, PCTreeRegistry<Key>> {
	PCTreeForest* m_pForest;

public:
	PCTreeRegistry(PCTreeForest* pcTreeForest) : m_pForest(pcTreeForest) { }

	//! Returns the index of \p key.
	static inline int keyToIndex(Key key);

	//! Returns whether \p key is associated with this registry.
	bool isKeyAssociated(Key key) const;

	//! Returns the maximum index of all keys managed by this registry.
	int maxKeyIndex() const;

	//! Returns the array size currently requested by this registry.
	int calculateArraySize(int add) const;

	operator PCTreeForest&() const { return *m_pForest; }

	operator PCTreeForest*() const { return m_pForest; }

	PCTreeForest& getForest() const { return *m_pForest; }
};

template<class Key>
bool PCTreeRegistry<Key>::isKeyAssociated(Key key) const {
#ifdef OGDF_DEBUG
	return key && key->getForest() == m_pForest;
#else
	return key;
#endif
}

template<class Key>
int PCTreeRegistry<Key>::keyToIndex(Key key) {
	return key->index();
}
}
