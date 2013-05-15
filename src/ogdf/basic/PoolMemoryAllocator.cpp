/*
 * $Revision: 3472 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-29 15:52:12 +0200 (Mo, 29. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of memory manager for more efficiently
 *        allocating small pieces of memory
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


#include <ogdf/basic/basic.h>


namespace ogdf {


struct PoolMemoryAllocator::PoolElement
{
	MemElemPtr m_gp;
	int        m_size;
};

struct PoolMemoryAllocator::BlockChain
{
	char m_fill[eBlockSize-sizeof(void*)];
	BlockChain *m_next;
};


PoolMemoryAllocator::PoolElement PoolMemoryAllocator::s_pool[eTableSize];
PoolMemoryAllocator::BlockChainPtr PoolMemoryAllocator::s_blocks;


#ifdef OGDF_MEMORY_POOL_NTS
PoolMemoryAllocator::MemElemPtr PoolMemoryAllocator::s_tp[eTableSize];

#elif defined(OGDF_NO_COMPILER_TLS)
CriticalSection *PoolMemoryAllocator::s_criticalSection;
pthread_key_t PoolMemoryAllocator::s_tpKey;

#else

CriticalSection *PoolMemoryAllocator::s_criticalSection;
OGDF_DECL_THREAD PoolMemoryAllocator::MemElemPtr PoolMemoryAllocator::s_tp[eTableSize];
#endif


void PoolMemoryAllocator::init()
{
#ifndef OGDF_MEMORY_POOL_NTS
#ifdef OGDF_NO_COMPILER_TLS
	pthread_key_create(&s_tpKey,NULL);
#endif

	s_criticalSection = new CriticalSection(500);
#endif

	initThread();
}


void PoolMemoryAllocator::initThread() {
#if !defined(OGDF_MEMORY_POOL_NTS) && defined(OGDF_NO_COMPILER_TLS)
		pthread_setspecific(s_tpKey,calloc(eTableSize,sizeof(MemElemPtr)));
#endif
}


void PoolMemoryAllocator::cleanup()
{
	BlockChainPtr p = s_blocks;
	while(p != 0) {
		BlockChainPtr pNext = p->m_next;
		free(p);
		p = pNext;
	}

#ifndef OGDF_MEMORY_POOL_NTS
#ifdef OGDF_NO_COMPILER_TLS
	pthread_key_delete(s_tpKey);
#endif

	delete s_criticalSection;
#endif
}


inline void PoolMemoryAllocator::enterCS()
{
#ifndef OGDF_MEMORY_POOL_NTS
	s_criticalSection->enter();
#endif
}


inline void PoolMemoryAllocator::leaveCS()
{
#ifndef OGDF_MEMORY_POOL_NTS
	s_criticalSection->leave();
#endif
}


bool PoolMemoryAllocator::checkSize(size_t nBytes) {
	return nBytes < eTableSize;
}


void *PoolMemoryAllocator::allocate(size_t nBytes) {
#if !defined(OGDF_MEMORY_POOL_NTS) && defined(OGDF_NO_COMPILER_TLS)
	MemElemPtr &pFreeBytes = *((MemElemPtr*)pthread_getspecific(s_tpKey)+nBytes);
#else
	MemElemPtr &pFreeBytes = s_tp[nBytes];
#endif
	if (OGDF_LIKELY(pFreeBytes != 0)) {
		MemElemPtr p = pFreeBytes;
		pFreeBytes = p->m_next;
		p->m_next = 0;
		return p;
	} else {
		return fillPool(pFreeBytes,__uint16(nBytes));
	}
}


void PoolMemoryAllocator::deallocate(size_t nBytes, void *p) {
#if !defined(OGDF_MEMORY_POOL_NTS) && defined(OGDF_NO_COMPILER_TLS)
	MemElemPtr &pFreeBytes = *((MemElemPtr*)pthread_getspecific(s_tpKey)+nBytes);
#else
	MemElemPtr &pFreeBytes = s_tp[nBytes];
#endif
	MemElemPtr(p)->m_next = pFreeBytes;
	pFreeBytes = MemElemPtr(p);
}


void PoolMemoryAllocator::deallocateList(size_t nBytes, void *pHead, void *pTail) {
#if !defined(OGDF_MEMORY_POOL_NTS) && defined(OGDF_NO_COMPILER_TLS)
	MemElemPtr &pFreeBytes = *((MemElemPtr*)pthread_getspecific(s_tpKey)+nBytes);
#else
	MemElemPtr &pFreeBytes = s_tp[nBytes];
#endif
	MemElemPtr(pTail)->m_next = pFreeBytes;
	pFreeBytes = MemElemPtr(pHead);
}



void PoolMemoryAllocator::flushPool()
{
#ifndef OGDF_MEMORY_POOL_NTS
	for(__uint16 nBytes = 1; nBytes < eTableSize; ++nBytes) {
#ifdef OGDF_NO_COMPILER_TLS
		MemElemPtr &pHead =  *((MemElemPtr*)pthread_getspecific(s_tpKey)+nBytes);
#else
		MemElemPtr &pHead = s_tp[nBytes];
#endif
		if(pHead != 0) {
			MemElemPtr pTail = pHead;
			int n = 1;

			while(pTail->m_next != 0) {
				pTail = pTail->m_next;
				++n;
			}

			MemElemPtr pOldHead = pHead;
			pHead = 0;

			enterCS();

			PoolElement &pe = s_pool[nBytes];

			pTail->m_next = pe.m_gp;
			pe.m_gp = pOldHead;
			pe.m_size += n;

			leaveCS();
		}
	}
#endif
}


void *PoolMemoryAllocator::fillPool(MemElemPtr &pFreeBytes, __uint16 nBytes)
{
	int nWords;
	int nSlices = slicesPerBlock(max(nBytes,(__uint16)eMinBytes),nWords);

#ifdef OGDF_MEMORY_POOL_NTS
	pFreeBytes = allocateBlock();
	makeSlices(pFreeBytes, nWords, nSlices);

#else
	enterCS();

	PoolElement &pe = s_pool[nBytes];
	if(pe.m_size >= nSlices) {
		MemElemPtr p = pFreeBytes = pe.m_gp;
		for(int i = 1; i < nSlices; ++i)
			p = p->m_next;

		pe.m_gp = p->m_next;
		pe.m_size -= nSlices;

		leaveCS();

		p->m_next = 0;

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


PoolMemoryAllocator::MemElemPtr
PoolMemoryAllocator::allocateBlock()
{
	BlockChainPtr pBlock = (BlockChainPtr) malloc(eBlockSize);

	pBlock->m_next = s_blocks;
	s_blocks = pBlock;

	return (MemElemPtr)pBlock;
}


void PoolMemoryAllocator::makeSlices(MemElemPtr pBlock, int nWords, int nSlices)
{
	do {
		pBlock = pBlock->m_next = pBlock+nWords;
	} while(--nSlices > 1);
	pBlock->m_next = 0;
}



size_t PoolMemoryAllocator::memoryAllocatedInBlocks()
{
	enterCS();

	size_t nBlocks = 0;
	for (BlockChainPtr p = s_blocks; p != 0; p = p->m_next)
		++nBlocks;

	leaveCS();

	return nBlocks * eBlockSize;
}


size_t PoolMemoryAllocator::memoryInGlobalFreeList()
{
	enterCS();

	size_t bytesFree = 0;
	for (int sz = 1; sz < eTableSize; ++sz)
	{
		const PoolElement &pe = s_pool[sz];
		bytesFree += pe.m_size * sz;
	}

	leaveCS();

	return bytesFree;
}


size_t PoolMemoryAllocator::memoryInThreadFreeList()
{
	size_t bytesFree = 0;
	for (int sz = 1; sz < eTableSize; ++sz)
	{
#if !defined(OGDF_MEMORY_POOL_NTS) && defined(OGDF_NO_COMPILER_TLS)
		MemElemPtr p = ((MemElemPtr*)pthread_getspecific(s_tpKey))[sz];
#else
		MemElemPtr p = s_tp[sz];
#endif
		for(; p != 0; p = p->m_next)
			bytesFree += sz;
	}

	return bytesFree;
}


void PoolMemoryAllocator::defrag()
{
	enterCS();

	int maxSize = 0;
	for(int sz = 1; sz < eTableSize; ++sz) {
		int size = s_pool[sz].m_size;
		maxSize = max(maxSize, size);
	}

	if(maxSize > 1) {
		MemElemPtr *a = new MemElemPtr[maxSize];

		for(int sz = 1; sz < eTableSize; ++sz)
		{
			PoolElement &pe = s_pool[sz];
			int n = pe.m_size;
			if(n > 1)
			{
				int i = 0;
				for(MemElemPtr p = pe.m_gp; p != 0; p = p->m_next)
					a[i++] = p;
				OGDF_ASSERT(i == n);
				std::sort(a, a+n);
				pe.m_gp = a[0];
				for(int i = 0; i < n-1; ++i) {
					a[i]->m_next = a[i+1];
				}
				a[n-1]->m_next = 0;
			}
		}
		delete [] a;
	}

	leaveCS();
}

}
