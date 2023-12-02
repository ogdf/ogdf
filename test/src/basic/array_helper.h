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
inline int maybeWrap(int value) {
	return value;
}

template<>
inline List<int> maybeWrap(int value) {
	return {value};
}

template<>
inline std::unique_ptr<int> maybeWrap(int value) {
	return std::unique_ptr<int>(new int(value));
}

template<>
inline std::vector<std::unique_ptr<int>> maybeWrap(int value) {
	std::vector<std::unique_ptr<int>> v;
	v.push_back(std::unique_ptr<int>(new int(value)));
	return v;
}

template<>
inline bool maybeWrap(int value) {
	return value;
}

template<typename Type>
inline int unwrap(Type& value);

template<>
inline int unwrap(std::unique_ptr<int>& value) {
	return *value;
}

template<>
inline int unwrap(std::vector<std::unique_ptr<int>>& value) {
	return *value.front();
}

/**
 * Perform basic tests for a map of registered keys to values.
 * @tparam BaseType the type of the class (e.g. Graph) that maintains the registered keys
 * @tparam ArrayType the type of array to be tested
 * @tparam KeyType the type of registered key
 * @tparam ElementType the value type
 *
 * @param title the title of the top-level bandit::describe
 * @param fillElement an arbitrary instance of \a ElementType
 * @param secondElement a second instance of \a ElementType, must differ from \p fillElement
 * @param initBase a function to initialize the base
 * @param chooseKey a function to choose an arbitrary key from the base
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key in the base
 */
template<class BaseType, template<typename, bool> class ArrayType, typename KeyType, typename ElementType>
void describeArray(const std::string& title, const ElementType& fillElement,
		const ElementType& secondElement, std::function<void(BaseType&)> initBase,
		std::function<KeyType(const BaseType&)> chooseKey,
		std::function<void(const BaseType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(BaseType&)> createKey) {
	using MyArrayType = ArrayType<ElementType, true>;
	using RegistryType = const typename MyArrayType::registry_type;
	using const_iterator = typename MyArrayType::const_iterator;
	using iterator = typename MyArrayType::iterator;

	describe(title.c_str(), [&]() {
		std::unique_ptr<MyArrayType> array;
		BaseType base;
		initBase(base);

		before_each([&]() {
			array.reset(new MyArrayType());
			initBase(base);
		});

		it("handles nested arrays well", [&]() {
			BaseType B;
			initBase(B);

			List<KeyType> keys;
			getAllKeys(B, keys);

			ArrayType<MyArrayType, true> nestedArray(B);
			for (KeyType k : keys) {
				nestedArray[k].init(B, fillElement);
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
				array->init(base);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
			});

			it("initializes w an empty registry", [&]() {
				BaseType B;
				array->init(B);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)B)));
				AssertThat(array->valid(), IsTrue());
			});

			it("initializes w a registry and filled", [&]() {
				array->init(base, fillElement);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(base)], Equals(fillElement));
			});

			it("is constructed w a registry", [&]() {
				array.reset(new MyArrayType(base));
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w an empty registry", [&]() {
				BaseType B;
				array.reset(new MyArrayType(B));
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)B)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w a registry and filled", [&]() {
				array.reset(new MyArrayType(base, fillElement));
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(base)], Equals(fillElement));
			});

			it("supports copy-construction", [&]() {
				array->init(base, fillElement);
				MyArrayType copiedArray(*array);
				AssertThat(copiedArray.registeredAt(), Equals(array->registeredAt()));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(base)], Equals(fillElement));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)], Equals(fillElement));
			});

			it("implements the assignment-operator", [&]() {
				array->init(base, fillElement);
				MyArrayType copiedArray = *array;
				AssertThat(copiedArray.registeredAt(), Equals(array->registeredAt()));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)], Equals(fillElement));
			});

			it("supports move-construction", [&]() {
				array->init(base, fillElement);
				MyArrayType copiedArray = std::move(*array);
				AssertThat(copiedArray.registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)], Equals(fillElement));
			});

			it("moves an array using the assignment operator", [&]() {
				array->init(base, fillElement);
				MyArrayType copiedArray;
				copiedArray = (std::move(*array));
				AssertThat(&(*copiedArray.registeredAt()), Equals(&((RegistryType&)base)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)], Equals(fillElement));
			});

			it("correctly fills the array with a default value", [&]() {
				array->init(base, fillElement);
				AssertThat(array->getDefault(), Equals(fillElement));
				array->setDefault(secondElement);
				array->fillWithDefault();
				AssertThat(array->getDefault(), Equals(secondElement));
				AssertThat((*array)[chooseKey(base)], Equals(secondElement));
			});

			it("allows filling the array with a specific element", [&]() {
				array->init(base, fillElement);
				array->fill(secondElement);
				AssertThat((*array)[chooseKey(base)], Equals(secondElement));
			});

			it("assigns the default value to a newly created key", [&]() {
				array->init(base, fillElement);
				KeyType key = createKey(base);
				AssertThat((*array)[key], Equals(fillElement));
			});

			it("retains the stored information when it grows", [&]() {
				array->init(base, fillElement);
				KeyType k = chooseKey(base);
				(*array)[k] = secondElement;
				int size = array->registeredAt()->getArraySize();

				for (int i = 0; i <= size; ++i) {
					createKey(base);
				}

				AssertThat(array->registeredAt()->getArraySize(), IsGreaterThan(size));
				AssertThat((*array)[k], Equals(secondElement));
			});
		});

		describe("access", [&]() {
			it("distinguishes between a valid and an invalid array", [&]() {
				AssertThat(array->valid(), IsFalse());
				array->init(base);
				AssertThat(array->valid(), IsTrue());
			});

			it("knows its registry", [&]() {
				array->init(base);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
			});

			it("allows access with the subscript operator", [&]() {
				array->init(base, fillElement);
				KeyType k = chooseKey(base);
				AssertThat((*array)[k], Equals(fillElement));
				AssertThat((*array)[k] = secondElement, Equals(secondElement));

				array->init(base, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray[k], Equals(fillElement));
			});

			it("allows access with the () operator", [&]() {
				array->init(base, fillElement);
				KeyType k = chooseKey(base);
				AssertThat((*array)(k), Equals(fillElement));
				AssertThat((*array)(k) = secondElement, Equals(secondElement));

				array->init(base, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray(k), Equals(fillElement));
			});
		});

		describe("iterators", [&]() {
			before_each([&]() { array->init(base, fillElement); });

			it("iterates over the array", [&]() {
				List<KeyType> list;
				getAllKeys(base, list);

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

/**
 * Perform basic tests for a map of registered keys to non copy constructible values.
 * @tparam BaseType the type of the class (e.g. Graph) that maintains the registered keys
 * @tparam ArrayType the type of array to be tested
 * @tparam KeyType the type of registered key
 * @tparam ElementType the value type
 *
 * @param title the title of the top-level bandit::describe
 * @param initBase a function to initialize the base
 * @param chooseKey a function to choose an arbitrary key from the base
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key in the base
 */
template<class BaseType, template<typename, bool> class ArrayType, typename KeyType, typename ElementType>
void describeArrayWithoutDefault(const std::string& title, std::function<void(BaseType&)> initBase,
		std::function<KeyType(const BaseType&)> chooseKey,
		std::function<void(const BaseType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(BaseType&)> createKey) {
	using MyArrayType = ArrayType<ElementType, false>;
	using RegistryType = const typename MyArrayType::registry_type;
	using const_iterator = typename MyArrayType::const_iterator;
	using iterator = typename MyArrayType::iterator;

	describe(title.c_str(), [&]() {
		std::unique_ptr<MyArrayType> array;
		BaseType base;
		initBase(base);

		before_each([&]() {
			array.reset(new MyArrayType());
			initBase(base);
		});

		it("handles nested arrays well", [&]() {
			BaseType B;
			initBase(B);

			List<KeyType> keys;
			getAllKeys(B, keys);

			ArrayType<MyArrayType, false> nestedArray(B);
			for (KeyType k : keys) {
				nestedArray[k].init(B);
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
				array->init(base);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
			});

			it("initializes w an empty registry", [&]() {
				BaseType B;
				array->init(B);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)B)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w a registry", [&]() {
				array.reset(new MyArrayType(base));
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
			});

			it("is constructed w an empty registry", [&]() {
				BaseType B;
				array.reset(new MyArrayType(B));
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)B)));
				AssertThat(array->valid(), IsTrue());
			});

			it("supports move-construction", [&]() {
				array->init(base);
				MyArrayType copiedArray = std::move(*array);
				AssertThat(copiedArray.registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)] == ElementType(), IsTrue());
			});

			it("moves an array using the assignment operator", [&]() {
				array->init(base);
				MyArrayType copiedArray;
				copiedArray = (std::move(*array));
				AssertThat(&(*copiedArray.registeredAt()), Equals(&((RegistryType&)base)));
				AssertThat(array->registeredAt(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(base)] == ElementType(), IsTrue());
			});

			it("assigns the default constructor value to a newly created key", [&]() {
				array->init(base);
				KeyType key = createKey(base);
				AssertThat((*array)[key] == ElementType(), IsTrue());
			});

			it("retains the stored information when it grows", [&]() {
				array->init(base);
				KeyType k = chooseKey(base);
				(*array)[k] = maybeWrap<ElementType>(42);
				int size = array->registeredAt()->getArraySize();

				for (int i = 0; i <= size; ++i) {
					createKey(base);
				}

				AssertThat(array->registeredAt()->getArraySize(), IsGreaterThan(size));
				AssertThat(unwrap<ElementType>((*array)[k]), Equals(42));
			});
		});

		describe("access", [&]() {
			it("distinguishes between a valid and an invalid array", [&]() {
				AssertThat(array->valid(), IsFalse());
				array->init(base);
				AssertThat(array->valid(), IsTrue());
			});

			it("knows its registry", [&]() {
				array->init(base);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
			});

			it("allows access with the subscript operator", [&]() {
				array->init(base);
				KeyType k = chooseKey(base);
				AssertThat((*array)[k] == ElementType(), IsTrue());
				(*array)[k] = maybeWrap<ElementType>(42);
				AssertThat(unwrap<ElementType>((*array)[k]), Equals(42));
			});

			it("allows access with the () operator", [&]() {
				array->init(base);
				KeyType k = chooseKey(base);
				AssertThat((*array)(k) == ElementType(), IsTrue());
				(*array)(k) = maybeWrap<ElementType>(42);
				AssertThat(unwrap<ElementType>((*array)(k)), Equals(42));
			});
		});

		describe("iterators", [&]() {
			before_each([&]() { array->init(base); });

			it("iterates over the array", [&]() {
				List<KeyType> list;
				getAllKeys(base, list);

				const MyArrayType cArray(base);
				int counter = 0;
				for (const_iterator it = cArray.begin(); it != cArray.end(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat((*array)[it.key()] == ElementType(), IsTrue());
					AssertThat(cArray[it.key()] == ElementType(), IsTrue());
					AssertThat(*it == ElementType(), IsTrue());
					(*array)[it.key()] = maybeWrap<ElementType>(counter);
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for (iterator it = array->begin(); it != array->end(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat(unwrap<ElementType>((*array)[it.key()]), Equals(counter));
					AssertThat(unwrap<ElementType>(*it), Equals(counter));
					*it = maybeWrap<ElementType>(counter + list.size());
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for (const_iterator it = cArray.cbegin(); it != cArray.cend(); it++) {
					AssertThat(counter, IsLessThan(list.size()));
					AssertThat(unwrap<ElementType>((*array)[it.key()]),
							Equals(counter + list.size()));
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				for (KeyType val : list) {
					AssertThat((*array)[val] != ElementType(), IsTrue());
				}
			});
		});
	});
}
