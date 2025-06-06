/*****************************************************************************************[Alloc.h]
Copyright (c) 2008-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#pragma once

#include <ogdf/basic/basic.h>
#include <ogdf/lib/minisat/mtl/XAlloc.h>
#include <ogdf/lib/minisat/mtl/Vec.h>

#pragma GCC visibility push(default)
namespace Minisat {
namespace Internal {

//=================================================================================================
// Simple Region-based memory allocator:

template<class T>
class RegionAllocator
{
	T*        memory;
	uint32_t  sz;
	uint32_t  cap;
	uint32_t  wasted_;

	void capacity(uint32_t min_cap);

 public:
	using Ref = uint32_t;

	enum { Ref_Undef = UINT_MAX };
	enum { Unit_Size = sizeof(uint32_t) };

	explicit RegionAllocator(uint32_t start_cap = 1024 * 1024) : memory(nullptr), sz(0), cap(0), wasted_(0){ capacity(start_cap); }

	~RegionAllocator()
	{
		if (memory != nullptr)
			::free(memory);
	}


	uint32_t size      () const      { return sz; }
	uint32_t wasted    () const      { return wasted_; }

	Ref      alloc     (int size);
	void     free      (int size)    { wasted_ += size; }

	// Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
	T&       operator[](Ref r)       { OGDF_ASSERT(r < sz); return memory[r]; }
	const T& operator[](Ref r) const { OGDF_ASSERT(r < sz); return memory[r]; }

	T*       lea       (Ref r)       { OGDF_ASSERT(r < sz); return &memory[r]; }
	const T* lea       (Ref r) const { OGDF_ASSERT(r < sz); return &memory[r]; }

	Ref ael(const T* t) {
		OGDF_ASSERT((const void*)t >= (const void*)&memory[0]);
		OGDF_ASSERT((const void*)t < (const void*)&memory[sz-1]);
		return (Ref)(t - &memory[0]);
	}

	void moveTo(RegionAllocator& to) {
		if (to.memory != nullptr) ::free(to.memory);
		to.memory = memory;
		to.sz = sz;
		to.cap = cap;
		to.wasted_ = wasted_;

		memory = nullptr;
		sz = cap = wasted_ = 0;
	}
};


template<class T>
void RegionAllocator<T>::capacity(uint32_t min_cap)
{
	if (cap >= min_cap) return;

	uint32_t prev_cap = cap;
	while (cap < min_cap) {
		// NOTE: Multiply by a factor (13/8) without causing overflow, then add 2 and make the
		// result even by clearing the least significant bit. The resulting sequence of capacities
		// is carefully chosen to hit a maximum capacity that is close to the '2^32-1' limit when
		// using 'uint32_t' as indices so that as much as possible of this space can be used.
		uint32_t delta = ((cap >> 1) + (cap >> 3) + 2) & ~1;
		cap += delta;

		if (cap <= prev_cap)
			throw OutOfMemoryException();
	}
#if 0
	 printf(" .. (%p) cap = %u\n", this, cap);
#endif

	OGDF_ASSERT(cap > 0);
	memory = static_cast<T*>(xrealloc(memory, sizeof(T)*cap));
}


template<class T>
typename RegionAllocator<T>::Ref
RegionAllocator<T>::alloc(int size)
{
#if 0
	 printf("ALLOC called (this = %p, size = %d)\n", this, size); fflush(stdout);
#endif
	OGDF_ASSERT(size > 0);
	capacity(sz + size);

	uint32_t prev_sz = sz;
	sz += size;

	// Handle overflow:
	if (sz < prev_sz)
		throw OutOfMemoryException();

	return prev_sz;
}


//=================================================================================================
} // namespace Internal
} // namespace Minisat
#pragma GCC visibility pop
