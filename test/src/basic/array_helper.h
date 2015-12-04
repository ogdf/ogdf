/** \file
 * \brief Generic tests for all array classes
 *
 * \author Stephan Beyer, Mirko Wagner, Tilo Wiedera
 *
 * \par License:
 * This file is part of the Open Graph Drawing, Drawing Framework (OGDF).
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
 * Software Foundation, Inc., 51 FRank,lin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/
#pragma once
#ifdef _MSC_VER
#pragma once
#endif


#include <bandit/bandit.h>
#include <ogdf/basic/graph_generators.h>

namespace ogdf
{

/**
 * Perform basic tests for a map of graph elements to values.
 * @tparam ArrayType the type of array to be tested
 * @tparam KeyType the type of key graph element
 * @tparam ElementType the value type
 *
 * @param title the title of the top-level bandit::describe
 * @param fillElement an arbitrary instance of \a ElementType
 * @param secondElement a second instance of \a ElementType, must differ from \a fillElement
 * @param chooseKey a function to choose an arbitrary key element from the graph
 * @param getAllKeys a function to generate a list of all keys
 * @param createKey a function to create a new key element in the graph
 */
template<template<typename> class ArrayType, typename KeyType, typename ElementType>
void describeArray(
		const std::string &title,
		const ElementType &fillElement,
		const ElementType &secondElement,
		std::function<KeyType (const Graph&)> chooseKey,
		std::function<void (const Graph&, List<KeyType>&)> getAllKeys,
		std::function<KeyType (Graph&)> createKey)
{
	using MyArrayType = ArrayType<ElementType>;
	using const_iterator = typename MyArrayType::const_iterator;
	using iterator = typename MyArrayType::iterator;

	bandit::describe(title.c_str(),[&]() {
		MyArrayType *array;
		Graph graph;
		graph = Graph();
		randomGraph(graph, 42, 168);

		bandit::before_each([&](){
			array = new MyArrayType();
		});

		bandit::after_each([&](){
			delete array;
		});

		bandit::it("handles nested arrays well", [&]() {
			Graph G;
			G.newEdge(G.newNode(), G.newNode());

			List<KeyType> keys;
			getAllKeys(G, keys);

			ArrayType<MyArrayType> nestedArray(G);
			for (KeyType k : keys) {
				nestedArray[k].init(G, fillElement);
			}
		});

		bandit::describe("init",[&](){
			bandit::it("initializes w/o a graph",[&](){
				AssertThat(array->graphOf(), IsNull());
				AssertThat(array->valid(), IsFalse());
				array->init();
				AssertThat(array->graphOf(), IsNull());
				AssertThat(array->valid(), IsFalse());
			});

			bandit::it("initializes w a graph",[&](){
				array->init(graph);
				AssertThat(array->graphOf(), Equals(&graph));
				AssertThat(array->valid(), IsTrue());
			});

			bandit::it("initializes w a graph and filled",[&](){
				array->init(graph, fillElement);
				AssertThat(array->graphOf(), Equals(&graph));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("is constructed w a graph",[&](){
				delete array;
				array = new MyArrayType(graph);
				AssertThat(array->graphOf(), Equals(&graph));
				AssertThat(array->valid(), IsTrue());
			});

			bandit::it("is constructed w a graph and filled",[&](){
				delete array;
				array = new MyArrayType(graph, fillElement);
				AssertThat(array->graphOf(), Equals(&graph));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("supports copy-construction",[&](){
				array->init(graph, fillElement);
				MyArrayType copiedArray(*array);
				AssertThat(copiedArray.graphOf(), Equals(array->graphOf()));
				AssertThat(array->valid(), IsTrue());
				AssertThat((*array)[chooseKey(graph)], Equals(fillElement));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("implements the assignment-operator",[&](){
				array->init(graph, fillElement);
				MyArrayType copiedArray = *array;
				AssertThat(copiedArray.graphOf(), Equals(array->graphOf()));
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("supports move-construction",[&](){
				array->init(graph, fillElement);
				MyArrayType copiedArray = std::move(*array);
				AssertThat(copiedArray.graphOf(), Equals(&graph));
				AssertThat(array->graphOf(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("moves an array using the assignment operator",[&](){
				array->init(graph, fillElement);
				MyArrayType copiedArray;
				copiedArray = (std::move(*array));
				AssertThat(&(*copiedArray.graphOf()), Equals(&graph));
				AssertThat(array->graphOf(), IsNull());
				AssertThat(array->valid(), IsFalse());
				AssertThat(copiedArray.valid(), IsTrue());
				AssertThat(copiedArray[chooseKey(graph)], Equals(fillElement));
			});

			bandit::it("assigns the default value to a newly created key", [&](){
				array->init(graph, fillElement);
				KeyType key = createKey(graph);
				AssertThat((*array)[key], Equals(fillElement));
			});
		});

		bandit::describe("access",[&](){
			bandit::it("distinguishes between a valid and an invalid array",[&](){
				AssertThat(array->valid(), IsFalse());
				array->init(graph);
				AssertThat(array->valid(), IsTrue());
			});

			bandit::it("knows its graph",[&](){
				array->init(graph);
				AssertThat(array->graphOf(), Equals(&graph));
			});

			bandit::it("allows access with the subscript operator",[&](){
				array->init(graph, fillElement);
				KeyType k = chooseKey(graph);
				AssertThat((*array)[k], Equals(fillElement));
				AssertThat((*array)[k] = secondElement, Equals(secondElement));

				array->init(graph, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray[k], Equals(fillElement));
			});

			bandit::it("allows access with the () operator",[&](){
				array->init(graph, fillElement);
				KeyType k = chooseKey(graph);
				AssertThat((*array)(k), Equals(fillElement));
				AssertThat((*array)(k) = secondElement, Equals(secondElement));

				array->init(graph, fillElement);
				const MyArrayType cAccessArray(*array);
				AssertThat(cAccessArray(k), Equals(fillElement));
			});
		});

		bandit::describe("iterators",[&](){
			bandit::before_each([&](){
				array->init(graph, fillElement);
			});

			bandit::it("iterates over the array",[&](){
				List<KeyType> list;
				getAllKeys(graph, list);

				const MyArrayType cArray(*array);
				int counter = 0;
				for(const_iterator it = cArray.begin(); it != cArray.end(); it++){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for(iterator it = array->begin(); it != array->end(); it++){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for(const_iterator it = cArray.cbegin(); it != cArray.cend(); it++){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));
			});

			bandit::it("iterates over the array backwards",[&](){
				List<KeyType> list;
				getAllKeys(graph, list);

				const MyArrayType cArray(*array);
				int counter = 0;
				for(const_iterator it = cArray.rbegin(); it != cArray.rend(); it--){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for(iterator it = array->rbegin(); it != array->rend(); it--){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));

				counter = 0;
				for(const_iterator it = cArray.crbegin(); it != cArray.crend(); it--){
					counter++;
				}
				AssertThat(counter, Equals(list.size()));
			});
		});
	});

}
}
