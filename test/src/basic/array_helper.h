/** \file
 * \brief Generic tests for all array classes
 *
 * \author Stephan Beyer, Mirko Wagner, Tilo Wiedera
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

#include <ogdf/basic/graph_generators.h>

#include <testing.h>

template<typename Type>
inline Type maybeWrap(int value);

template<>
inline int maybeWrap<int>(int value) {
	return value;
}

template<>
inline List<int> maybeWrap<List<int>>(int value) {
	return {value};
}

/**
 * Perform basic tests for a map of registry elements to values.
 * @tparam RegistryType the type of the associated registry
 * @tparam ArrayType the type of array to be tested
 * @tparam KeyType the type of key registry element
 * @tparam ElementType the value type
 *
 * @param title the title of the top-level bandit::describe
 * @param fillElement an arbitrary instance of \a ElementType
 * @param secondElement a second instance of \a ElementType, must differ from \p fillElement
 * @param initRegistry a function to initialize the registry
 * @param chooseKey a function to choose an arbitrary key element from the registry
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key element in the registry
 */
template<class RegistryType, template<typename> class ArrayType, typename KeyType, typename ElementType>
void describeArray(const std::string& title, const ElementType& fillElement,
		const ElementType& secondElement, std::function<void(RegistryType&)> initRegistry,
		std::function<KeyType(const RegistryType&)> chooseKey,
		std::function<void(const RegistryType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(RegistryType&)> createKey) {
	using MyArrayType = ArrayType<ElementType>;
	using RegistryBaseType = const typename MyArrayType::registry_type;
	using const_iterator = typename MyArrayType::const_iterator;
	using iterator = typename MyArrayType::iterator;

	describe(title.c_str(), [&]() {
		std::unique_ptr<MyArrayType> array;
		RegistryType registry;
		initRegistry(registry);

		before_each([&]() {
			array.reset(new MyArrayType());
			initRegistry(registry);
		});

		it("handles nested arrays well", [&]() {
			RegistryType R;
			initRegistry(R);

			List<KeyType> keys;
			getAllKeys(R, keys);

			ArrayType<MyArrayType> nestedArray(R);
			for (KeyType k : keys) {
				nestedArray[k].init(R, fillElement);
			}
		});

		describe("init", [&]() {
			it("initializes w/o a registry", [&]() {
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				array->init();
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
			});

			it("initializes w a registry", [&]() {
				array->init(registry);
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->valid(), IsTrue());
			});

			it("initializes w an empty registry", [&]() {
				RegistryType R;
				array->init(R);
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)R)));
				AssertThat(array->valid(), IsTrue());
			});

			it("initializes w a registry and filled", [&]() {
				array->init(registry, fillElement);
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(registry)], Equals(fillElement));
			});

			it("is constructed w a registry", [&]() {
				array.reset(new MyArrayType(registry));
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w an empty registry", [&]() {
				RegistryType R;
				array.reset(new MyArrayType(R));
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)R)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w a registry and filled", [&]() {
				array.reset(new MyArrayType(registry, fillElement));
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(registry)], Equals(fillElement));
			});

			it("supports copy-construction", [&]() {
				array->init(registry, fillElement);
				MyArrayType copiedArray(*array);
				AssertThat(copiedArray.registeredAt(), Equals(array->registeredAt()));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(registry)], Equals(fillElement));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(registry)], Equals(fillElement));
			});

			it("implements the assignment-operator", [&]() {
				array->init(registry, fillElement);
				MyArrayType copiedArray = *array;
				AssertThat(copiedArray.registeredAt(), Equals(array->registeredAt()));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(registry)], Equals(fillElement));
			});

			it("supports move-construction", [&]() {
				array->init(registry, fillElement);
				MyArrayType copiedArray = std::move(*array);
				AssertThat(copiedArray.registeredAt(), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(registry)], Equals(fillElement));
			});

			it("moves an array using the assignment operator", [&]() {
				array->init(registry, fillElement);
				MyArrayType copiedArray;
				copiedArray = (std::move(*array));
				AssertThat(&(*copiedArray.registeredAt()), Equals(&((RegistryBaseType&)registry)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(registry)], Equals(fillElement));
			});

			it("correctly fills the array with a default value", [&]() {
				array->init(registry, fillElement);
				AssertThat(array->getDefault(), Equals(fillElement));
				array->setDefault(secondElement);
				array->fillWithDefault();
				AssertThat(array->getDefault(), Equals(secondElement));
				AssertThat((*array)[chooseKey(registry)], Equals(secondElement));
			});

			it("allows filling the array with a specific element", [&]() {
				array->init(registry, fillElement);
				array->fill(secondElement);
				AssertThat((*array)[chooseKey(registry)], Equals(secondElement));
			});

			it("assigns the default value to a newly created key", [&]() {
				array->init(registry, fillElement);
				KeyType key = createKey(registry);
				AssertThat((*array)[key], Equals(fillElement));
			});

			it("retains the stored information when it grows", [&]() {
				array->init(registry, fillElement);
				KeyType k = chooseKey(registry);
				(*array)[k] = secondElement;
				int size = array->registeredAt()->getArraySize();

				for (int i = 0; i <= size; ++i) {
					createKey(registry);
				}

				AssertThat(array->registeredAt()->getArraySize(), IsGreaterThan(size));
				AssertThat((*array)[k], Equals(secondElement));
			});
		});

		describe("access", [&]() {
			it("distinguishes between a valid and an invalid array", [&]() {
				AssertThat(array->valid(), IsFalse());
				array->init(registry);
				AssertThat(array->valid(), IsTrue());
			});

			it("knows its registry", [&]() {
				array->init(registry);
				AssertThat(array->registeredAt(), Equals(&((RegistryBaseType&)registry)));
			});

			it("allows access with the subscript operator", [&]() {
				array->init(registry, fillElement);
				KeyType k = chooseKey(registry);
				AssertThat((*array)[k], Equals(fillElement));
				AssertThat((*array)[k] = secondElement, Equals(secondElement));

				array->init(registry, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray[k], Equals(fillElement));
			});

			it("allows access with the () operator", [&]() {
				array->init(registry, fillElement);
				KeyType k = chooseKey(registry);
				AssertThat((*array)(k), Equals(fillElement));
				AssertThat((*array)(k) = secondElement, Equals(secondElement));

				array->init(registry, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray(k), Equals(fillElement));
			});
		});

		describe("iterators", [&]() {
			before_each([&]() { array->init(registry, fillElement); });

			it("iterates over the array", [&]() {
				List<KeyType> list;
				getAllKeys(registry, list);

				const MyArrayType cArray(*array);
				int counter = 0;
				for (const_iterator it = cArray.begin(); it != cArray.end(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat((*array)[it.key()], Equals(fillElement));
					AssertThat(cArray[it.key()], Equals(fillElement));
					AssertThat(*it, Equals(fillElement));
					(*array)[it.key()] = maybeWrap<ElementType>(counter);
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				int defaults = 0;
				for (iterator it = array->begin(); it != array->end(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat((*array)[it.key()], Equals(maybeWrap<ElementType>(counter)));
					AssertThat(*it, Equals(maybeWrap<ElementType>(counter)));
					*it = maybeWrap<ElementType>(counter + list.size());
					if (maybeWrap<ElementType>(counter + list.size()) == fillElement) {
						defaults++;
					}
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for (const_iterator it = cArray.cbegin(); it != cArray.cend(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat((*array)[it.key()],
							Equals(maybeWrap<ElementType>(counter + list.size())));
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				for (KeyType val : list) {
					if ((*array)[val] == fillElement) {
						defaults--;
					}
				}
				AssertThat(defaults, Equals(0));
			});
		});
	});
}
