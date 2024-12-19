/** \file
 * \brief Tests for ogdf::PoolMemoryAllocator, ogdf::MallocMemoryAllocator and the respective macros.
 *
 * \author Tilo Wiedera
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

#include <ogdf/basic/System.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/memory.h>
#include <ogdf/basic/memory/PoolMemoryAllocator.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include <testing.h>

template<size_t size>
class OGDFObject {
	char x[size];

	OGDF_NEW_DELETE
};

template<size_t size>
class MallocObject {
	char x[size];

	OGDF_MALLOC_NEW_DELETE
};

template<template<size_t> class ObjectOfSize>
void describeMemoryManager(const std::string& name) {
	describe(name + " allocator", [] {
		after_each([] {
			size_t managerAllocated = System::memoryAllocatedByMemoryManager();
			size_t mallocAllocated = System::memoryAllocatedByMalloc();
			AssertThat(managerAllocated, IsLessThanOrEqualTo(mallocAllocated));
		});

		it("allocates objects that need exactly 1 byte", [] { delete new ObjectOfSize<1>; });

		it("does not deallocate nullptr", [] {
			ObjectOfSize<1>* ptr = nullptr;
			delete ptr;
		});
	});
}

go_bandit([] {
	describeMemoryManager<OGDFObject>("OGDF");
	describeMemoryManager<MallocObject>("Malloc");

#define OBJ_SIZE 128

	describe(
			"PoolMemoryManager",
			[] {
				std::vector<OGDFObject<OBJ_SIZE>*> objects;
				const size_t MAX = 1024;
				const size_t ALL_OBJS_SIZE = MAX * sizeof(OGDFObject<OBJ_SIZE>);
				objects.reserve(MAX);

				it("allocates memory", [&] {
					auto get_used = [] {
						return OGDF_ALLOCATOR::memoryAllocatedInBlocks()
								- OGDF_ALLOCATOR::memoryInGlobalFreeList()
								- OGDF_ALLOCATOR::memoryInThreadFreeList();
					};
					size_t used_before = get_used();
					for (int i = 0; i < MAX; ++i) {
						objects.push_back(new OGDFObject<OBJ_SIZE>());
					}
					AssertThat(get_used() - used_before, IsGreaterThanOrEqualTo(ALL_OBJS_SIZE));
				});

				it("grows its free list when objects are deleted", [&] {
					std::vector<size_t> sizes_old;
					PoolMemoryAllocator::getThreadFreeListSizes(sizes_old);

					size_t mAllocatedInBlocks = PoolMemoryAllocator::memoryAllocatedInBlocks();
					size_t mInGlobalFreeList = PoolMemoryAllocator::memoryInGlobalFreeList();
					size_t mInThreadFreeList = PoolMemoryAllocator::memoryInThreadFreeList();
					std::shuffle(objects.begin(), objects.end(), std::mt19937(randomSeed()));
					for (auto object : objects) {
						delete object;
					}
					AssertThat(PoolMemoryAllocator::memoryAllocatedInBlocks(),
							Equals(mAllocatedInBlocks));
					AssertThat(PoolMemoryAllocator::memoryInGlobalFreeList(),
							Equals(mInGlobalFreeList));
					AssertThat(PoolMemoryAllocator::memoryInThreadFreeList(),
							IsGreaterThan(mInThreadFreeList));
					AssertThat(PoolMemoryAllocator::memoryInThreadFreeList() - mInThreadFreeList,
							Equals(ALL_OBJS_SIZE));

					std::vector<size_t> sizes_new;
					PoolMemoryAllocator::getThreadFreeListSizes(sizes_new);
					AssertThat(sizes_new.at(OBJ_SIZE) - sizes_old.at(OBJ_SIZE), Equals(MAX));
					sizes_new[OBJ_SIZE] = sizes_old[OBJ_SIZE];
					AssertThat(sizes_new, Equals(sizes_old));

					objects.clear();
				});

				it("correctly defragments its thread free list", [&] {
					PoolMemoryAllocator::defragThread();
					for (int i = 0; i < MAX; ++i) { // first set of objects from main thread for later
						objects.push_back(new OGDFObject<OBJ_SIZE>());
					}
					AssertThat(std::is_sorted(objects.begin(), objects.end()), IsTrue());
				});

				std::vector<size_t> thread_sizes_old;
				PoolMemoryAllocator::getThreadFreeListSizes(thread_sizes_old);
				std::vector<size_t> global_sizes_old;
				PoolMemoryAllocator::getGlobalFreeListSizes(global_sizes_old);
				std::vector<size_t> side_thread_sizes_new;

				ogdf::Thread([&] {
					it("handles allocations and deletions mixed between threads", [&] {
						for (int i = 0; i < MAX; ++i) { // second set of objects from side thread
							objects.push_back(new OGDFObject<OBJ_SIZE>());
						}
						std::shuffle(objects.begin(), objects.end(), std::mt19937(randomSeed()));
						AssertThat(objects.size(), Equals(MAX * 2));
						for (int i = 0; i < MAX; ++i) { // half of objects deleted in side thread
							delete objects.back();
							objects.pop_back();
						}
						AssertThat(objects.size(), Equals(MAX));
						PoolMemoryAllocator::getThreadFreeListSizes(side_thread_sizes_new);
					});
				}).join();
				it("handles deletions of objects from other threads", [&] {
					for (auto object : objects) { // last half of objects deleted in main thread
						delete object;
					}
					objects.clear();
				});

				it("hands thread memory over to global pool on thread termination", [&] {
					std::vector<size_t> thread_sizes_new;
					PoolMemoryAllocator::getThreadFreeListSizes(thread_sizes_new);
					std::vector<size_t> global_sizes_new;
					PoolMemoryAllocator::getGlobalFreeListSizes(global_sizes_new);

#ifdef OGDF_MEMORY_POOL_NTS
					AssertThat(thread_sizes_new.at(OBJ_SIZE), IsGreaterThanOrEqualTo(MAX * 2));
					AssertThat(global_sizes_new.empty(), IsTrue());
#else
			AssertThat(thread_sizes_new.at(OBJ_SIZE),
					Equals(thread_sizes_old.at(OBJ_SIZE) + MAX));
			AssertThat(global_sizes_new.at(OBJ_SIZE),
					Equals(global_sizes_old.at(OBJ_SIZE) + side_thread_sizes_new.at(OBJ_SIZE)));
#endif
				});


				it(
						"correctly defragments its global free list",
						[&] {
							// update counters
							PoolMemoryAllocator::getGlobalFreeListSizes(global_sizes_old);
							size_t mAllocatedInBlocks =
									PoolMemoryAllocator::memoryAllocatedInBlocks();

							PoolMemoryAllocator::defragGlobal(); // thread list will not be defrag'ed...

							// so use up those guys first
							std::vector<size_t> thread_sizes_new;
							PoolMemoryAllocator::getThreadFreeListSizes(thread_sizes_new);
							std::vector<OGDFObject<OBJ_SIZE>*> objects2;
							for (int i = 0; i < thread_sizes_new.at(OBJ_SIZE); ++i) {
								objects2.push_back(new OGDFObject<OBJ_SIZE>());
							}
							PoolMemoryAllocator::getThreadFreeListSizes(thread_sizes_new);
							AssertThat(thread_sizes_new.at(OBJ_SIZE), Equals(0));

							// now we should be using from the global list
							for (int i = 0; i < MAX; ++i) {
								objects.push_back(new OGDFObject<OBJ_SIZE>());
							}
							AssertThat(std::is_sorted(objects.begin(), objects.end()), IsTrue());

							std::vector<size_t> global_sizes_new;
							PoolMemoryAllocator::getGlobalFreeListSizes(global_sizes_new);
							AssertThat(global_sizes_new.at(OBJ_SIZE) - global_sizes_old.at(OBJ_SIZE),
									IsGreaterThanOrEqualTo(MAX));

							AssertThat(PoolMemoryAllocator::memoryAllocatedInBlocks(),
									Equals(mAllocatedInBlocks));

							for (auto object : objects) {
								delete object;
							}
							for (auto object : objects2) {
								delete object;
							}
						},
#ifdef OGDF_MEMORY_POOL_NTS
						true
#else
						false
#endif
				);
			},
#ifdef OGDF_MEMORY_MALLOC_TS
			true
#else
			false
#endif
	);
});
