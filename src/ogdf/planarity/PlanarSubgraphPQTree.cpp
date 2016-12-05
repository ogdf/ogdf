/** \file
 * \brief Implementation of the class PlanarSubgraphPQTree.
 *
 * Implements a PQTree with added features for the planarity test.
 * Used by BoothLueker.
 *
 * \author Sebastian Leipert
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

#include <ogdf/internal/planarity/PlanarSubgraphPQTree.h>

namespace ogdf{

// Replaces the pertinent subtree by a P-node with leaves as children
// corresponding to the incoming edges of the node v. These edges
// are to be specified by their keys stored in leafKeys.
void PlanarSubgraphPQTree::
ReplaceRoot(SListPure<PlanarLeafKey<whaInfo*>*> &leafKeys)
{
	if (m_pertinentRoot->status() == PQNodeRoot::FULL)
		ReplaceFullRoot(leafKeys);
	else
		ReplacePartialRoot(leafKeys);
}

// Initializes a PQTree by a set of leaves that will korrespond to
// the set of Keys stored in leafKeys.
int PlanarSubgraphPQTree::
Initialize(SListPure<PlanarLeafKey<whaInfo*>*> &leafKeys)
{
	SListPure<PQLeafKey<edge, whaInfo*, bool>*> castLeafKeys;
	for (PlanarLeafKey<whaInfo*> *leafPtr : leafKeys)
		castLeafKeys.pushBack(static_cast<PQLeafKey<edge, whaInfo*, bool>*>(leafPtr));

	return PQTree<edge, whaInfo*, bool>::Initialize(castLeafKeys);
}


// Reduction reduced a set of leaves determined by their keys stored
// in leafKeys. Integer redNumber is for debugging only.
bool PlanarSubgraphPQTree::Reduction(
	SListPure<PlanarLeafKey<whaInfo*>*>     &leafKeys,
	SList<PQLeafKey<edge, whaInfo*, bool>*> &eliminatedKeys)
{
	SListPure<PQLeafKey<edge, whaInfo*, bool>*> castLeafKeys;

	for (PlanarLeafKey<whaInfo*>*leafPtr : leafKeys)
	{
		castLeafKeys.pushBack(static_cast<PQLeafKey<edge, whaInfo*, bool>*>(leafPtr));
	}

	determineMinRemoveSequence(castLeafKeys, eliminatedKeys);
	removeEliminatedLeaves(eliminatedKeys);

	SListIterator<PQLeafKey<edge,whaInfo*,bool>* >  itn = castLeafKeys.begin();
	SListIterator<PQLeafKey<edge,whaInfo*,bool>* >  itp = itn++;
	for (; itn.valid();)
	{
		if ((*itn)->nodePointer()->status()== PQNodeRoot::WHA_DELETE)
		{
			++itn;
			castLeafKeys.delSucc(itp);
		}
		else
			itp = itn++;
	}

	if ((*castLeafKeys.begin())->nodePointer()->status() == PQNodeRoot::WHA_DELETE)
		castLeafKeys.popFront();


	return Reduce(castLeafKeys);
}



// Function ReplaceFullRoot either replaces the full root
// or one full child of a partial root of a pertinent subtree
// by a single P-node  with leaves corresponding the keys stored in leafKeys.
void PlanarSubgraphPQTree::
ReplaceFullRoot(SListPure<PlanarLeafKey<whaInfo*>*> &leafKeys)
{

	PQLeaf<edge, whaInfo*, bool>          *leafPtr = nullptr; // dummy
	PQInternalNode<edge, whaInfo*, bool>	*nodePtr = nullptr; // dummy
	PQNode<edge, whaInfo*, bool>		    *currentNode = nullptr; // dummy

	if (!leafKeys.empty() && leafKeys.front() == leafKeys.back())
	{
		//ReplaceFullRoot: replace pertinent root by a single leaf
		leafPtr = new PQLeaf<edge, whaInfo*, bool>(m_identificationNumber++,
			PQNodeRoot::EMPTY, (PQLeafKey<edge, whaInfo*, bool>*)leafKeys.front());
		exchangeNodes(m_pertinentRoot, (PQNode<edge, whaInfo*, bool>*) leafPtr);
		if (m_pertinentRoot == m_root)
			m_root = (PQNode<edge, whaInfo*, bool>*) leafPtr;
	}
	else if (!leafKeys.empty()) // at least two leaves
	{
		//replace pertinent root by a $P$-node
		if ((m_pertinentRoot->type() == PQNodeRoot::PNode) ||
			(m_pertinentRoot->type() == PQNodeRoot::QNode))
		{
			nodePtr = (PQInternalNode<edge, whaInfo*, bool>*)m_pertinentRoot;
			nodePtr->type(PQNodeRoot::PNode);
			nodePtr->status(PQNodeRoot::PERTROOT);
			nodePtr->childCount(0);
			while (!fullChildren(m_pertinentRoot)->empty())
			{
				currentNode = fullChildren(m_pertinentRoot)->popFrontRet();
				removeChildFromSiblings(currentNode);
			}
		}
		else if (m_pertinentRoot->type() == PQNodeRoot::leaf)
		{
			nodePtr = new PQInternalNode<edge, whaInfo*, bool>(m_identificationNumber++,
				PQNodeRoot::PNode, PQNodeRoot::EMPTY);
			exchangeNodes(m_pertinentRoot, nodePtr);
		}

		SListPure<PQLeafKey<edge, whaInfo*, bool>*> castLeafKeys;
		for (PlanarLeafKey<whaInfo*>* leafPtr : leafKeys)
			castLeafKeys.pushBack(static_cast<PQLeafKey<edge, whaInfo*, bool>*>(leafPtr));
		addNewLeavesToTree(nodePtr, castLeafKeys);
	}

}


// Function ReplacePartialRoot replaces all full nodes by a single P-node
// with leaves corresponding the keys stored in leafKeys.
void PlanarSubgraphPQTree::
	ReplacePartialRoot(SListPure<PlanarLeafKey<whaInfo*>*> &leafKeys)

{
	PQNode<edge,whaInfo*,bool>  *currentNode = nullptr;

	m_pertinentRoot->childCount(m_pertinentRoot->childCount() + 1 -
		fullChildren(m_pertinentRoot)->size());

	while (fullChildren(m_pertinentRoot)->size() > 1)
	{
		currentNode = fullChildren(m_pertinentRoot)->popFrontRet();
		removeChildFromSiblings(currentNode);
	}

	currentNode = fullChildren(m_pertinentRoot)->popFrontRet();

	currentNode->parent(m_pertinentRoot);
	m_pertinentRoot = currentNode;
	ReplaceFullRoot(leafKeys);

}


/**
The function removeEliminatedLeaves handles the difficult task of
cleaning up after every reduction.

After a reduction is complete, different kind of garbage has to be handled.
- Pertinent leaves that are not in the maximal pertinent sequence.
  from the $PQ$-tree in order to get it reducable have to be deleted.
- The memory of some pertinent nodes, that have only pertinent leaves not beeing
  in the maximal pertinent sequence in their frontier has to be freed.
- Pertinent nodes that have only one child left after the removal
  of pertinent leaves not beeing in the maximal pertinent sequence
  have to be deleted.
- The memory of all full nodes has to be freed, since the complete
  pertinent subtree is replaced by a $P$-node after the reduction.
- Nodes, that have been removed during the call of the function [[Reduce]]
  of the base class template [[PQTree]] from the $PQ$-tree have to be
  kept but marked as nonexisting.
*/

/**************************************************************************************
							removeEliminatedLeaves
***************************************************************************************/

void PlanarSubgraphPQTree::
removeEliminatedLeaves(SList<PQLeafKey<edge, whaInfo*, bool>*> &eliminatedKeys)
{
	for (PQLeafKey<edge, whaInfo*, bool> *key : eliminatedKeys)
	{
		PQNode<edge, whaInfo*, bool>* nodePtr = key->nodePointer();
		PQNode<edge, whaInfo*, bool>* parent = nodePtr->parent();
		PQNode<edge, whaInfo*, bool>* sibling = nodePtr->getNextSib(nullptr);

		removeNodeFromTree(parent, nodePtr);
		checkIfOnlyChild(sibling, parent);
		if (parent->status() == PQNodeRoot::TO_BE_DELETED)
		{
			parent->status(PQNodeRoot::WHA_DELETE);
		}
		nodePtr->status(PQNodeRoot::WHA_DELETE);
	}
}


}
