/** \file
 * \brief Implementation of memory manager for more efficiently
 *        allocating small pieces of memory
 *
 * \author Carsten Gutwenger, Tilo Wiedera
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory/MallocMemoryAllocator.h>
#include <ogdf/basic/memory/PoolMemoryAllocator.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <vector>

namespace ogdf {

struct PoolMemoryAllocator::PoolElement {
	MemElemPtr m_gp;
	int m_size;
};

struct PoolMemoryAllocator::BlockChain {
	char m_fill[BLOCK_SIZE - sizeof(void*)];
	BlockChain* m_next;
};

PoolMemoryAllocator::BlockChain* PoolMemoryAllocator::s_blocks;

#ifdef OGDF_DEBUG
long long PoolMemoryAllocator::s_globallyAllocatedBytes = 0;
thread_local long long PoolMemoryAllocator::s_locallyAllocatedBytes = 0;
#endif

#ifdef OGDF_MEMORY_POOL_NTS
PoolMemoryAllocator::MemElemPtr PoolMemoryAllocator::s_tp[TABLE_SIZE];
#else
PoolMemoryAllocator::PoolElement PoolMemoryAllocator::s_pool[TABLE_SIZE];
std::mutex PoolMemoryAllocator::s_mutex;
thread_local PoolMemoryAllocator::MemElemPtr PoolMemoryAllocator::s_tp[TABLE_SIZE];
#endif

void PoolMemoryAllocator::cleanup() {
	// check if all memory is correctly freed (if not we have a memory leak)
	OGDF_ASSERT(s_globallyAllocatedBytes + s_locallyAllocatedBytes == 0);

	BlockChain* p = s_blocks;
	while (p != nullptr) {
		BlockChain* pNext = p->m_next;
		free(p);
		p = pNext;
	}
}

bool PoolMemoryAllocator::checkSize(size_t nBytes) { return nBytes < TABLE_SIZE; }

void* PoolMemoryAllocator::allocate(size_t nBytes) {
	MemElemPtr& pFreeBytes = s_tp[nBytes];
	void* result;

	if (OGDF_LIKELY(pFreeBytes != nullptr)) {
		MemElemPtr p = pFreeBytes;
		pFreeBytes = p->m_next;
		p->m_next = nullptr;
		result = p;
	} else {
		result = fillPool(pFreeBytes, uint16_t(nBytes));
	}
#ifdef OGDF_DEBUG
	s_locallyAllocatedBytes += (result == nullptr ? 0 : nBytes);
#endif
	return result;
}

void PoolMemoryAllocator::deallocate(size_t nBytes, void* p) {
	MemElemPtr& pFreeBytes = s_tp[nBytes];
	MemElemPtr(p)->m_next = pFreeBytes;
	pFreeBytes = MemElemPtr(p);
#ifdef OGDF_DEBUG
	s_locallyAllocatedBytes -= nBytes;
#endif
}

void PoolMemoryAllocator::deallocateList(size_t nBytes, void* pHead, void* pTail) {
	if (OGDF_LIKELY(nBytes < TABLE_SIZE)) {
#ifdef OGDF_DEBUG
		for (auto ptr = MemElemPtr(pHead); ptr != MemElemPtr(pTail)->m_next; ptr = ptr->m_next) {
			s_locallyAllocatedBytes -= nBytes;
		}
#endif
		MemElemPtr& pFreeBytes = s_tp[nBytes];
		MemElemPtr(pTail)->m_next = pFreeBytes;
		pFreeBytes = MemElemPtr(pHead);
	} else {
		MallocMemoryAllocator::deallocateList(nBytes, pHead, pTail);
	}
}

void PoolMemoryAllocator::flushPool() {
#ifndef OGDF_MEMORY_POOL_NTS
	for (uint16_t nBytes = 1; nBytes < TABLE_SIZE; ++nBytes) {
		MemElemPtr& pHead = s_tp[nBytes];
		if (pHead != nullptr) {
			MemElemPtr pTail = pHead;
			int n = 1;

			while (pTail->m_next != nullptr) {
				pTail = pTail->m_next;
				++n;
			}

			MemElemPtr pOldHead = pHead;
			pHead = nullptr;

			enterCS();

			PoolElement& pe = s_pool[nBytes];

			pTail->m_next = pe.m_gp;
			pe.m_gp = pOldHead;
			pe.m_size += n;

			leaveCS();
		}
	}
#	ifdef OGDF_DEBUG
	enterCS();
	s_globallyAllocatedBytes += s_locallyAllocatedBytes;
	leaveCS();
	s_locallyAllocatedBytes = 0;
#	endif
#endif
}

void* PoolMemoryAllocator::fillPool(MemElemPtr& pFreeBytes, uint16_t nBytes) {
	int nWords;
	int nSlices = slicesPerBlock(max(nBytes, (uint16_t)MIN_BYTES), nWords);

#ifdef OGDF_MEMORY_POOL_NTS
	pFreeBytes = allocateBlock();
	makeSlices(pFreeBytes, nWords, nSlices);
#else
	enterCS();

	PoolElement& pe = s_pool[nBytes];
	// only take elements from global pool if it contains at least as many
	// as we would get from a new block
	if (pe.m_size >= nSlices) {
		MemElemPtr p = pFreeBytes = pe.m_gp;
		for (int i = 1; i < nSlices; ++i) {
			p = p->m_next;
		}

		pe.m_gp = p->m_next;
		pe.m_size -= nSlices;

		leaveCS();

		p->m_next = nullptr;

	} else {
		pFreeBytes = allocateBlock();

		leaveCS();

		makeSlices(pFreeBytes, nWords, nSlices);
	}
#endif

	MemElemPtr p = pFreeBytes;
	pFreeBytes = p->m_next;
	return p;
}

PoolMemoryAllocator::MemElemPtr PoolMemoryAllocator::allocateBlock() {
	BlockChain* pBlock = static_cast<BlockChain*>(malloc(BLOCK_SIZE));

	pBlock->m_next = s_blocks;
	s_blocks = pBlock;

	return reinterpret_cast<MemElemPtr>(pBlock);
}

void PoolMemoryAllocator::makeSlices(MemElemPtr pBlock, int nWords, int nSlices) {
	do {
		pBlock = pBlock->m_next = pBlock + nWords;
	} while (--nSlices > 1);
	pBlock->m_next = nullptr;
}

size_t PoolMemoryAllocator::memoryAllocatedInBlocks() {
	enterCS();

	size_t nBlocks = 0;
	for (BlockChain* p = s_blocks; p != nullptr; p = p->m_next) {
		++nBlocks;
	}

	leaveCS();

	return nBlocks * BLOCK_SIZE;
}

size_t PoolMemoryAllocator::memoryInGlobalFreeList() {
	size_t bytesFree = 0;
#ifndef OGDF_MEMORY_POOL_NTS
	enterCS();
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		const PoolElement& pe = s_pool[sz];
		bytesFree += pe.m_size * sz;
	}
	leaveCS();
#endif

	return bytesFree;
}

size_t PoolMemoryAllocator::memoryInThreadFreeList() {
	size_t bytesFree = 0;
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		for (MemElemPtr p = s_tp[sz]; p != nullptr; p = p->m_next) {
			bytesFree += sz;
		}
	}
	return bytesFree;
}

void PoolMemoryAllocator::defragGlobal() {
#ifndef OGDF_MEMORY_POOL_NTS
	enterCS();

	std::vector<MemElemPtr> elems;
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		PoolElement& pe = s_pool[sz];
		elems.reserve(pe.m_size);
		for (auto p = pe.m_gp; p != nullptr; p = p->m_next) {
			elems.push_back(p);
		}
		OGDF_ASSERT(elems.size() == pe.m_size);
		if (elems.empty()) {
			continue;
		}
		std::sort(elems.begin(), elems.end());

		MemElemPtr pred = pe.m_gp = elems.front();
		for (auto p : elems) {
			pred->m_next = p;
			pred = p;
		}
		elems.back()->m_next = nullptr;
		elems.clear();
	}

	leaveCS();
#endif
}

void PoolMemoryAllocator::defragThread() {
	std::vector<MemElemPtr> elems;
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		for (auto p = s_tp[sz]; p != nullptr; p = p->m_next) {
			elems.push_back(p);
		}
		if (elems.empty()) {
			continue;
		}
		std::sort(elems.begin(), elems.end());

		MemElemPtr pred = s_tp[sz] = elems.front();
		for (auto p : elems) {
			pred->m_next = p;
			pred = p;
		}
		elems.back()->m_next = nullptr;
		elems.clear();
	}
}

void PoolMemoryAllocator::getGlobalFreeListSizes(std::vector<size_t>& sizes) {
#ifdef OGDF_MEMORY_POOL_NTS
	sizes.clear();
#else
	enterCS();
	sizes.reserve(TABLE_SIZE);
	sizes.assign(1, 0);
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		sizes.push_back(s_pool[sz].m_size);
	}
	leaveCS();
#endif
}

void PoolMemoryAllocator::getThreadFreeListSizes(std::vector<size_t>& sizes) {
	sizes.reserve(TABLE_SIZE);
	sizes.assign(1, 0);
	for (size_t sz = 1; sz < TABLE_SIZE; ++sz) {
		size_t size = 0;
		for (auto p = s_tp[sz]; p != nullptr; p = p->m_next) {
			size++;
		}
		sizes.push_back(size);
	}
}

}
