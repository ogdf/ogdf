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

#include <ogdf/basic/List.h>

#include <functional>
#include <memory>
#include <string>
#include <vector>

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

template<class BaseType, class MyArrayType, class RegistryType>
void runInitTests(const BaseType& base, std::unique_ptr<MyArrayType>& array) {
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
			runInitTests<BaseType, MyArrayType, RegistryType>(base, array);

			it("initializes w a registry and filled", [&]() {
				array->init(base, fillElement);
				AssertThat(array->registeredAt(), Equals(&((RegistryType&)base)));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(base)], Equals(fillElement));
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
			runInitTests<BaseType, MyArrayType, RegistryType>(base, array);

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

/**
 * Perform several tests for a type of registered array.
 *
 * @tparam BaseType the type of the class (e.g. Graph) that maintains the registered keys
 * @tparam ArrayType the type of array to be tested
 * @tparam KeyType the type of registered key
 *
 * @param arrayType the name of the array type
 * @param initBase a function to initialize the base
 * @param chooseKey a function to choose an arbitrary key from the base
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key in the base
 */
template<class BaseType, template<typename, bool> class ArrayType, typename KeyType>
void runBasicArrayTests(const std::string& arrayType, std::function<void(BaseType&)> initBase,
		std::function<KeyType(const BaseType&)> chooseKey,
		std::function<void(const BaseType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(BaseType&)> createKey) {
	describeArray<BaseType, ArrayType, KeyType, int>( //
			arrayType + " filled with ints", //
			42, 43, //
			initBase, chooseKey, getAllKeys, createKey);
	describeArray<BaseType, ArrayType, KeyType, List<int>>( //
			arrayType + " filled with lists of ints", //
			{1, 2, 3}, {42}, //
			initBase, chooseKey, getAllKeys, createKey);
	describeArray<BaseType, ArrayType, KeyType, bool>( //
			arrayType + " filled with bools", //
			false, true, //
			initBase, chooseKey, getAllKeys, createKey);

	describeArrayWithoutDefault<BaseType, ArrayType, KeyType, std::unique_ptr<int>>( //
			arrayType + " filled with unique pointers", //
			initBase, chooseKey, getAllKeys, createKey);
	describeArrayWithoutDefault<BaseType, ArrayType, KeyType, std::vector<std::unique_ptr<int>>>( //
			arrayType + " filled with vectors of unique pointers", //
			initBase, chooseKey, getAllKeys, createKey);
}

template<class BaseType, typename SetType, typename KeyType>
void describeSet(const std::string& setType, std::function<void(BaseType&)> initBase,
		std::function<KeyType(const BaseType&)> chooseKey,
		std::function<void(const BaseType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(BaseType&)> createKey,
		std::function<void(BaseType&, KeyType)> deleteKey,
		std::function<void(BaseType&)> clearAllKeys, std::function<bool(KeyType, KeyType)> equals) {
	describe(setType, [&]() {
		it("initializes with an empty Registry", [&]() {
			BaseType B;
			SetType S(B);
			AssertThat(S.size(), Equals(0));
		});

		{
			BaseType B;
			initBase(B);
			it("initializes with a non-empty Registry", [&]() {
				SetType S(B);
				AssertThat(S.size(), Equals(0));
			});

			KeyType a = chooseKey(B);
			KeyType b = chooseKey(B);
			while (equals(a, b)) {
				b = chooseKey(B);
			}
			it("allows adding values", [&]() {
				SetType S(B);
				AssertThat(S.isMember(a), IsFalse());
				AssertThat(S.isMember(b), IsFalse());
				S.insert(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.elements().size(), Equals(1));
				AssertThat(S.elements().front(), Equals(a));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsFalse());
			});
			it("allows adding values twice", [&]() {
				SetType S(B);
				S.insert(a);
				S.insert(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsFalse());
			});
			it("allows removing and adding values back", [&]() {
				SetType S(B);
				S.insert(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.isMember(a), IsTrue());
				S.remove(a);
				AssertThat(S.size(), Equals(0));
				AssertThat(S.isMember(a), IsFalse());

				S.remove(a);
				AssertThat(S.size(), Equals(0));
				AssertThat(S.isMember(a), IsFalse());

				S.insert(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsFalse());
			});
			it("allows clearing", [&]() {
				SetType S(B);
				S.insert(a);
				S.insert(b);
				AssertThat(S.size(), Equals(2));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsTrue());

				S.clear();
				AssertThat(S.size(), Equals(0));
				AssertThat(S.elements().empty(), IsTrue());
				AssertThat(S.isMember(a), IsFalse());

				S.remove(a);
				AssertThat(S.size(), Equals(0));
				AssertThat(S.isMember(a), IsFalse());

				S.insert(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsFalse());
			});
			SetType onlyA(B);
			onlyA.insert(a);
			it("allows copying", [&]() {
				SetType S(onlyA);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.elements().front(), Equals(a));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsFalse());

				S.insert(b);
				AssertThat(S.size(), Equals(2));
				AssertThat(S.isMember(a), IsTrue());
				AssertThat(S.isMember(b), IsTrue());

				S.remove(a);
				AssertThat(S.size(), Equals(1));
				AssertThat(S.elements().front(), Equals(b));
				AssertThat(S.isMember(a), IsFalse());
				AssertThat(S.isMember(b), IsTrue());
			});
			it("checks equality", [&]() {
				SetType S1(onlyA);
				SetType S2(B);
				S2.insert(a);

				AssertThat(S1, Equals(S2));
				AssertThat(S1, Equals(onlyA));
				AssertThat(S2, Equals(onlyA));
				AssertThat(onlyA, Equals(S1));
				AssertThat(onlyA, Equals(S2));

				S1.insert(a);
				AssertThat(S1, Equals(S2));
				AssertThat(S1, Equals(onlyA));
				AssertThat(S2, Equals(onlyA));
				AssertThat(onlyA, Equals(S1));
				AssertThat(onlyA, Equals(S2));

				S1.insert(b);
				AssertThat(S1, !Equals(S2));
				AssertThat(S1, !Equals(onlyA));
				AssertThat(S2, Equals(onlyA));
				AssertThat(onlyA, !Equals(S1));
				AssertThat(onlyA, Equals(S2));
			});
		}

		it("is notified of key deletion", [&]() {
			BaseType B;
			initBase(B);

			KeyType a = chooseKey(B);
			KeyType b = chooseKey(B);
			while (equals(a, b)) {
				b = chooseKey(B);
			}

			SetType S(B);
			S.insert(a);
			S.insert(b);
			AssertThat(S.size(), Equals(2));
			deleteKey(B, a);
			AssertThat(S.size(), Equals(1));
			AssertThat(S.elements().front(), Equals(b));
		});
		it("is notified of key clearing", [&]() {
			BaseType B;
			initBase(B);
			SetType S(B);

			List<KeyType> list;
			getAllKeys(B, list);
			for (KeyType k : list) {
				AssertThat(S.isMember(k), IsFalse());
				S.insert(k);
				AssertThat(S.isMember(k), IsTrue());
			}
			AssertThat(S.size(), Equals(list.size()));
			S.clear();
			AssertThat(S.size(), Equals(0));

			for (KeyType k : list) {
				AssertThat(S.isMember(k), IsFalse());
				S.insert(k);
				AssertThat(S.isMember(k), IsTrue());
			}
			AssertThat(S.size(), Equals(list.size()));

			clearAllKeys(B);
			AssertThat(S.size(), Equals(0));
			AssertThat(S.elements().empty(), IsTrue());
		});
		it("allows changing its base", [&]() {
			std::unique_ptr<BaseType> A = std::make_unique<BaseType>();
			std::unique_ptr<BaseType> B = std::make_unique<BaseType>();
			initBase(*A);
			initBase(*B);
			SetType S(*A);

#define OGDF_CAST_REGISTRY_PTR(ptr) (&static_cast<const typename SetType::registry_type&>(*(ptr)))

			KeyType a = chooseKey(*A);
			S.insert(a);
			AssertThat(S.size(), Equals(1));
			AssertThat(OGDF_CAST_REGISTRY_PTR(S.registeredAt()), Equals(OGDF_CAST_REGISTRY_PTR(A)));
			S.init(*B);
			AssertThat(OGDF_CAST_REGISTRY_PTR(S.registeredAt()), Equals(OGDF_CAST_REGISTRY_PTR(B)));
			AssertThat(S.size(), Equals(0));
			KeyType b = chooseKey(*B);
			S.insert(b);
			AssertThat(S.size(), Equals(1));

			deleteKey(*A, a);
			AssertThat(S.size(), Equals(1));
			deleteKey(*B, b);
			AssertThat(S.size(), Equals(0));
			KeyType b2 = chooseKey(*B);
			S.insert(b2);
			AssertThat(S.size(), Equals(1));

			A.reset();
			AssertThat(S.size(), Equals(1));
			AssertThat(S.elements().front(), Equals(b2));

			B.reset();
			AssertThat(S.size(), Equals(0));
			AssertThat(S.registeredAt(), IsNull());
		});

		it("allows equality comparisons", [&]() {
			std::unique_ptr<BaseType> G = std::make_unique<BaseType>();
			initBase(*G);

			SetType SA(*G);
			SetType SB;
			AssertThat(SA, !Equals(SB));
			SB.init(*G);
			AssertThat(SA, Equals(SB));

			KeyType n = chooseKey(*G);
			SA.insert(n);
			AssertThat(SA, !Equals(SB));
			SB.insert(n);
			AssertThat(SA, Equals(SB));

			KeyType a = n;
			while (equals(a, n)) {
				a = chooseKey(*G);
			}
			KeyType b = n;
			while (equals(n, b) || equals(a, b)) {
				b = chooseKey(*G);
			}

			SA.insert(a);
			AssertThat(SA, !Equals(SB));
			SB.insert(b);
			AssertThat(SA, !Equals(SB));
			SB.insert(a);
			AssertThat(SA, !Equals(SB));
			SA.insert(b);
			AssertThat(SA, Equals(SB));

			SA.remove(n);
			AssertThat(SA, !Equals(SB));
			SB.remove(n);
			AssertThat(SA, Equals(SB));

			SA.init();
			AssertThat(SA, !Equals(SB));
			G.reset();
			AssertThat(SA, Equals(SB));
		});
		it("allows copy construction and assignment", [&]() {
			std::unique_ptr<BaseType> G = std::make_unique<BaseType>();
			initBase(*G);

			KeyType a = chooseKey(*G);
			KeyType b = a;
			while (equals(a, b)) {
				b = chooseKey(*G);
			}

			SetType SA(*G);
			SA.insert(a);
			SA.insert(b);
			SetType SB(SA);
			AssertThat(SB.isMember(a), IsTrue());
			AssertThat(SB.isMember(b), IsTrue());
			AssertThat(SB.size(), Equals(2));
			AssertThat(SA, Equals(SB));
			AssertThat(SB, Equals(SA));

			SB.remove(a);
			AssertThat(SA.isMember(a), IsTrue());
			AssertThat(SA.isMember(b), IsTrue());
			AssertThat(SA.size(), Equals(2));
			AssertThat(SB.isMember(a), IsFalse());
			AssertThat(SB.isMember(b), IsTrue());
			AssertThat(SB.size(), Equals(1));
			AssertThat(SA != SB, IsTrue());
			AssertThat(SB != SA, IsTrue());

			SB = SA;
			AssertThat(SA.isMember(a), IsTrue());
			AssertThat(SA.isMember(b), IsTrue());
			AssertThat(SA.size(), Equals(2));
			AssertThat(SB.isMember(a), IsTrue());
			AssertThat(SB.isMember(b), IsTrue());
			AssertThat(SB.size(), Equals(2));
			AssertThat(SA, Equals(SB));
			AssertThat(SB, Equals(SA));

			SB = SetType(*G);
			AssertThat(SA.isMember(a), IsTrue());
			AssertThat(SA.isMember(b), IsTrue());
			AssertThat(SA.size(), Equals(2));
			AssertThat(SB.isMember(a), IsFalse());
			AssertThat(SB.isMember(b), IsFalse());
			AssertThat(SB.size(), Equals(0));
			AssertThat(SA != SB, IsTrue());
			AssertThat(SB != SA, IsTrue());
		});
	});
};

/**
 * Perform several tests for a type of registered set.
 *
 * @tparam BaseType the type of the class (e.g. Graph) that maintains the registered keys
 * @tparam SetType the type of set to be tested
 * @tparam KeyType the type of registered key
 *
 * @param setType the name of the set type
 * @param initBase a function to initialize the base
 * @param chooseKey a function to choose an arbitrary key from the base
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key in the base
 * @param deleteKey a function to delete an existing key in the base
 * @param clearAllKeys a function to delete all keys in the base
 * @param equals a function to check whether two objects are equal wrt. being deleted together
 */
template<class BaseType, class SetType, typename KeyType>
void runBasicSetTests(
		const std::string& setType, std::function<void(BaseType&)> initBase,
		std::function<KeyType(const BaseType&)> chooseKey,
		std::function<void(const BaseType&, List<KeyType>&)> getAllKeys,
		std::function<KeyType(BaseType&)> createKey,
		std::function<void(BaseType&, KeyType)> deleteKey,
		std::function<void(BaseType&)> clearAllKeys,
		std::function<bool(KeyType, KeyType)> equals = [](KeyType a, KeyType b) { return a == b; }) {
	describeSet<BaseType, SetType, KeyType>(setType, initBase, chooseKey, getAllKeys, createKey,
			deleteKey, clearAllKeys, equals);
	// FIXME is cleared() called before or after the changes?
}
