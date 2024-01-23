#pragma once

#include "PCEnum.h"

namespace pc_tree {
class PCTreeForest;

template<class Key>
class PCTreeRegistry : public ogdf::RegistryBase<Key, PCTreeRegistry<Key>> {
	PCTreeForest* m_pForest;

public:
	PCTreeRegistry(PCTreeForest* pcTreeForest) : m_pForest(pcTreeForest) { }

	bool isKeyAssociated(Key key) const override;

	int keyToIndex(Key key) const override;

	int keyArrayTableSize() const override;

	int maxKeyIndex() const override;

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
int PCTreeRegistry<Key>::keyToIndex(Key key) const {
	return key->index();
}
}
