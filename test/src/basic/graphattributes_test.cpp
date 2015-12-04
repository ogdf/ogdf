/** \file
 * \brief test of class GraphAttributes.
 *
 * Class GraphAttributes extends a graph by graphical attributes like
 * node position, color, etc.
 *
 * \author Mirko Wagner
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
#include <bandit/bandit.h>

#include <ogdf/basic/DualGraph.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <resources.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/graphics.h>

using namespace ogdf;
using namespace bandit;

using Type = ogdf::GraphAttributes;
GraphAttributes grattr;

/**
 * tests getter and setter of an attribute
 *
 * \param functionToTest is both a setter and getter
 * \param constFunctionToTest is a const getter
 * \param defaultValue is the value the attribute is initialized to
 * \param secondValue is the value the attribute gets set to
 * \param neededAttributes the attributes which are needed for the attribute
 * \param attributeName a string of the name
 * \param_t testElementType the type which the attribute is
 */
template<typename testElementType>
void testAttribute(
	std::function<testElementType& (void)> functionToTest,
	std::function<testElementType (void)> constFunctionToTest,
	testElementType defaultValue,
	testElementType secondValue,
	long neededAttributes,
	string attributeName)
{
	describe(attributeName,[&](){
		before_each([&](){
			grattr.initAttributes(neededAttributes);
		});

		it("throws an Exception, if "+attributeName+" is called without the right GraphAttributes set",[&](){
#ifdef OGDF_DEBUG
			grattr.destroyAttributes(neededAttributes);
			AssertThrows(PreconditionViolatedException, functionToTest());
#endif
		});

		it("gets "+attributeName,[&](){
			AssertThat(constFunctionToTest(), Equals(defaultValue));
			AssertThat(functionToTest(), Equals(defaultValue));
		});

		it("sets "+attributeName,[&](){
			testElementType& value = functionToTest();
			value = secondValue;
			AssertThat(functionToTest(), Equals(secondValue));
			AssertThat(constFunctionToTest(), Equals(secondValue));
		});
	});
}

go_bandit([](){
describe("graph attributes",[&](){
	Graph graph;
	long attr =0;

	before_each([&](){
		graph = Graph();
	});

	it("initializes w/o a graph",[&](){
		grattr = GraphAttributes();
		AssertThat(&grattr.constGraph(), IsNull());
		AssertThat(grattr.attributes(), Equals(0));
	});

	it("initializes w a graph and w attributes and the constructor",[&](){
		attr = Type::nodeId;
		grattr = GraphAttributes(graph, attr);
		AssertThat(&grattr.constGraph(), Equals(&graph));
		AssertThat(grattr.attributes(), Equals(Type::nodeId));
	});

	it("initializes w a graph and w/o attributes and the constructor",[&](){
		grattr = GraphAttributes(graph);
		AssertThat(&grattr.constGraph(), Equals(&graph));
		AssertThat(grattr.attributes(), Equals(Type::nodeGraphics | Type::edgeGraphics));
	});

	it("initializes w a graph and init",[&](){
		grattr = GraphAttributes();
		attr = Type::nodeId;
		grattr.init(graph, attr);
		AssertThat(&grattr.constGraph(), Equals(&graph));
		AssertThat(grattr.attributes(), Equals(Type::nodeId));
	});

	it("destroys its attributes",[&](){
		attr = Type::nodeGraphics | Type::nodeLabel;
		grattr = GraphAttributes(graph, attr);
		AssertThat(&grattr.constGraph(), Equals(&graph));
		AssertThat(grattr.attributes(), Equals(Type::nodeGraphics | Type::nodeLabel));
		grattr.destroyAttributes(Type::nodeGraphics | Type::nodeId);
		AssertThat(grattr.attributes(), Equals(Type::nodeLabel));
	});

	it("initializes its attributes",[&](){
		attr = Type::nodeGraphics | Type::nodeLabel;
		grattr = GraphAttributes(graph, attr);
		AssertThat(&grattr.constGraph(), Equals(&graph));
		AssertThat(grattr.attributes(), Equals(Type::nodeGraphics | Type::nodeLabel));
		grattr.initAttributes(Type::nodeId | Type::nodeLabel);
		AssertThat(grattr.attributes(), Equals(Type::nodeGraphics | Type::nodeLabel | Type::nodeId));
	});

	it("knows its current Attributes",[&](){
		attr = Type::nodeId | Type::nodeLabel;
		grattr = GraphAttributes(graph, attr);
		AssertThat(grattr.has(attr), IsTrue());
		AssertThat(grattr.has(Type::nodeId), IsTrue());
		AssertThat(grattr.has(Type::nodeId | Type::nodeGraphics), IsFalse());
		AssertThat(grattr.has(Type::nodeGraphics), IsFalse());
	});

	describe("attributes",[&](){
		node v;
		edge e;
		const GraphAttributes &cGrattr = grattr;

		before_each([&](){
			graph = Graph();
			completeGraph(graph, 7);
			grattr = GraphAttributes(graph, 0);
			v = graph.chooseNode();
			e = graph.chooseEdge();
		});

		it("knows if it's directed",[&](){
			AssertThat(grattr.directed(), IsTrue());
			grattr.setDirected(false);
			AssertThat(grattr.directed(), IsFalse());
		});

		testAttribute<double>(
			[&]() -> double& { return grattr.x(v); },
			[&](){ return cGrattr.x(v); },
			0, 42,
			Type::nodeGraphics, "x");

		testAttribute<double>(
			[&]() -> double& { return grattr.xLabel(v); },
			[&](){ return grattr.xLabel(v); },
			0, 42,
			Type::nodeLabel | Type::nodeLabelPosition, "xLabel");

		testAttribute<double>(
			[&]() -> double& { return grattr.y(v); },
			[&](){ return cGrattr.y(v); },
			0, 42,
			Type::nodeGraphics, "y");

		testAttribute<double>(
			[&]() -> double& { return grattr.yLabel(v); },
			[&](){ return cGrattr.yLabel(v); },
			0, 42,
			Type::nodeLabel | Type::nodeLabelPosition, "yLabel");

		testAttribute<double>(
			[&]() -> double& { return grattr.z(v); },
			[&](){ return cGrattr.z(v); },
			0, 42,
			Type::nodeGraphics | Type::threeD, "z");

		testAttribute<double>(
			[&]() -> double& { return grattr.zLabel(v); },
			[&](){ return cGrattr.zLabel(v); },
			0, 42,
			Type::nodeLabel | Type::nodeLabelPosition | Type::threeD | Type::nodeGraphics, "zLabel");

		testAttribute<double>(
			[&]() -> double& { return grattr.width(v); },
			[&](){ return cGrattr.width(v); },
			20, 42,
			Type::nodeGraphics, "width of a node");

		it("it assings width using a NodeArray",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.width());
#endif
			grattr.initAttributes(Type::nodeGraphics);
			NodeArray<double> widthNA(graph,42);
			AssertThat(cGrattr.width().graphOf(), Equals(&graph));
			AssertThat(grattr.width().graphOf(), Equals(&graph));
			AssertThat(cGrattr.width()[v], Equals(LayoutStandards::defaultNodeWidth()));
			AssertThat(grattr.width()[v], Equals(LayoutStandards::defaultNodeWidth()));
			grattr.width() = widthNA;
			AssertThat(cGrattr.width()[v], Equals(42));
			AssertThat(grattr.width()[v], Equals(42));
			grattr.setAllWidth(1337);
			AssertThat(cGrattr.width(v), Equals(1337));
			AssertThat(grattr.width(v), Equals(1337));
		});

		testAttribute<double>(
			[&]() -> double& { return grattr.height(v); },
			[&](){ return cGrattr.height(v); },
			20, 42,
			Type::nodeGraphics, "height of a node");

		it("it assigns height using a NodeArray",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.height());
#endif
			grattr.initAttributes(Type::nodeGraphics);
			NodeArray<double> heightNA(graph,42);
			AssertThat(cGrattr.height().graphOf(), Equals(&graph));
			AssertThat(grattr.height().graphOf(), Equals(&graph));
			grattr.height() = heightNA;
			AssertThat(cGrattr.height()[v], Equals(42));
			AssertThat(grattr.height()[v], Equals(42));
			grattr.setAllHeight(1337);
			AssertThat(cGrattr.height(v), Equals(1337));
			AssertThat(grattr.height(v), Equals(1337));
		});

		testAttribute<int>(
			[&]() -> int& { return grattr.weight(v); },
			[&](){ return cGrattr.weight(v); },
			0, 42,
			Type::nodeWeight, "weight of a node");

		testAttribute<Graph::EdgeType>(
			[&]() -> Graph::EdgeType& { return grattr.type(e); },
			[&](){ return cGrattr.type(e); },
			Graph::association, Graph::generalization,
			Type::edgeType, "type of an edge");

		testAttribute<Graph::NodeType>(
			[&]() -> Graph::NodeType& { return grattr.type(v); },
			[&](){ return cGrattr.type(v); },
			Graph::vertex, Graph::dummy,
			Type::nodeType, "type of a node");

		testAttribute<uint32_t>(
			[&]() -> uint32_t& { return grattr.subGraphBits(e); },
			[&](){ return cGrattr.subGraphBits(e); },
			0, 42,
			Type::edgeSubGraphs, "SubGraphBits");

		testAttribute<float>(
			[&]() -> float& { return grattr.strokeWidth(e); },
			[&](){ return cGrattr.strokeWidth(e); },
			LayoutStandards::defaultEdgeStroke().m_width, 42,
			Type::edgeStyle | Type::edgeGraphics, "strokeWidth edge");

		testAttribute<float>(
			[&]() -> float& { return grattr.strokeWidth(v); },
			[&](){ return cGrattr.strokeWidth(v); },
			LayoutStandards::defaultEdgeStroke().m_width, 42,
			Type::nodeStyle | Type::nodeGraphics, "strokeWidth node");

		it("strokeType node",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.strokeType(v));
#endif
			grattr.initAttributes(Type::nodeStyle | Type::nodeGraphics);
			AssertThat(cGrattr.strokeType(v), Equals(LayoutStandards::defaultNodeStroke().m_type));
			AssertThat(grattr.strokeType(v), Equals(LayoutStandards::defaultNodeStroke().m_type));
			grattr.setStrokeType(v, StrokeType::stDot);
			AssertThat(grattr.strokeType(v), Equals(StrokeType::stDot));
		});

		it("strokeType edge",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.strokeType(e));
#endif
			grattr.initAttributes(Type::edgeStyle | Type::edgeGraphics);
			AssertThat(cGrattr.strokeType(e), Equals(LayoutStandards::defaultEdgeStroke().m_type));
			AssertThat(grattr.strokeType(e), Equals(LayoutStandards::defaultEdgeStroke().m_type));
			grattr.setStrokeType(e, StrokeType::stDot);
			AssertThat(grattr.strokeType(e), Equals(StrokeType::stDot));
		});

		testAttribute<Color>(
			[&]() -> Color& { return grattr.strokeColor(e); },
			[&](){ return cGrattr.strokeColor(e); },
			LayoutStandards::defaultEdgeStroke().m_color, ogdf::Color::Name::Turquoise,
			Type::edgeStyle | Type::edgeGraphics, "strokeColor edge");

		testAttribute<Color>(
			[&]() -> Color& { return grattr.strokeColor(v); },
			[&](){ return cGrattr.strokeColor(v); },
			LayoutStandards::defaultEdgeStroke().m_color, ogdf::Color::Name::Turquoise,
			Type::nodeStyle | Type::nodeGraphics, "strokeColor node");

		testAttribute<Shape>(
			[&]() -> Shape& { return grattr.shape(v); },
			[&](){ return cGrattr.shape(v); },
			LayoutStandards::defaultNodeShape(), Shape::shRect,
			Type::nodeGraphics, "shape node");

		testAttribute<EdgeArrow>(
			[&]() -> EdgeArrow& { return grattr.arrowType(e); },
			[&](){ return cGrattr.arrowType(e); },
			LayoutStandards::defaultEdgeArrow(), EdgeArrow::eaBoth,
			Type::edgeArrow, "arrowType");

		testAttribute<double>(
			[&]() -> double& { return grattr.doubleWeight(e); },
			[&](){ return cGrattr.doubleWeight(e); },
			1.0, 42.0,
			Type::edgeDoubleWeight, "doubleWeight");

		testAttribute<Color>(
			[&]() -> Color& { return grattr.fillBgColor(v); },
			[&](){ return cGrattr.fillBgColor(v); },
			LayoutStandards::defaultNodeFill().m_bgColor, ogdf::Color::Turquoise,
			Type::nodeStyle | Type::nodeGraphics, "fillBgColor");

		testAttribute<Color>(
			[&]() -> Color& { return grattr.fillColor(v); },
			[&](){ return cGrattr.fillColor(v); },
			LayoutStandards::defaultNodeFill().m_color,  Color(ogdf::Color::Name::Turquoise),
			Type::nodeStyle | Type::nodeGraphics, "fillColor");

		it("fillPattern",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.fillPattern(v));
#endif
			grattr.initAttributes(Type::nodeStyle | Type::nodeGraphics);
			AssertThat(cGrattr.fillPattern(v), Equals(LayoutStandards::defaultNodeFill().m_pattern));
			AssertThat(grattr.fillPattern(v), Equals(LayoutStandards::defaultNodeFill().m_pattern));
			grattr.setFillPattern(v, FillPattern::fpCross);
			AssertThat(grattr.fillPattern(v), Equals(FillPattern::fpCross));
		});

		grattr.initAttributes(Type::nodeId);
		testAttribute<int>(
			[&]() -> int& { return grattr.idNode(v); },
			[&]() { return cGrattr.idNode(v); },
			grattr.idNode(v), 42,
			Type::nodeId, "idNode");
		grattr.destroyAttributes(Type::nodeId);

		it("(in|add|remove)SubGraph",[&](){
#ifdef OGDF_DEBUG
			AssertThrows(PreconditionViolatedException, grattr.inSubGraph(e, 13));
#endif
			grattr.initAttributes(Type::edgeSubGraphs);
			AssertThat(cGrattr.inSubGraph(e, 13), IsFalse());
			AssertThat(grattr.inSubGraph(e, 13), IsFalse());
			grattr.addSubGraph(e, 13);
			AssertThat(cGrattr.inSubGraph(e, 13), IsTrue());
			AssertThat(grattr.inSubGraph(e, 13), IsTrue());
			grattr.removeSubGraph(e, 13);
			AssertThat(cGrattr.inSubGraph(e, 13), IsFalse());
			AssertThat(grattr.inSubGraph(e, 13), IsFalse());
		});

		testAttribute<int>(
			[&]() -> int& { return grattr.intWeight(e); },
			[&](){ return cGrattr.intWeight(e); },
			1,  42,
			Type::edgeIntWeight, "intWeight");

		testAttribute<string>(
			[&]() -> string& { return grattr.label(v); },
			[&](){ return cGrattr.label(v); },
			"",  "ogdf" ,
			Type::nodeLabel, "label node");

		testAttribute<string>(
			[&]() -> string& { return grattr.label(e); },
			[&](){ return cGrattr.label(e); },
			"",  "ogdf" ,
			Type::edgeLabel, "label edge");


		testAttribute<string>(
			[&]() -> string& { return grattr.templateNode(v); },
			[&](){ return cGrattr.templateNode(v); },
			"",  "ogdf" ,
			Type::nodeTemplate, "template node");
	});

	describe("change position of elements",[&](){
		before_each([&](){
			graph = Graph();
			completeGraph(graph, 100);
			grattr = GraphAttributes(graph, 0);
			grattr.initAttributes(Type::nodeGraphics | Type::edgeGraphics);
			for(node v : graph.nodes){
				grattr.x(v) = ogdf::randomNumber(-100, 100);
				grattr.y(v) = ogdf::randomNumber(-100, 100);
			}
			grattr.addNodeCenter2Bends(1);
			grattr.translateToNonNeg();
		});

		after_each([&](){
		});

		it("translates to non-negative coordinates",[&](){
			for(node v : graph.nodes){
				AssertThat(grattr.x(v) - grattr.width(v) / 2, IsGreaterThan(0) || Equals(0));
				AssertThat(grattr.y(v) - grattr.width(v) / 2, IsGreaterThan(0) || Equals(0));
			}
			for(edge e : graph.edges){
				for(DPoint &p : grattr.bends(e)){
					AssertThat(p.m_x, IsGreaterThan(0) || Equals(0));
					AssertThat(p.m_y, IsGreaterThan(0) || Equals(0));
				}
			}
		});

		it("translates",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.translate(1.0, 42.0);
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(ga.x(v)+1.0));
				AssertThat(grattr.y(v), Equals(ga.y(v)+42.0));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(p_old.m_x+1.0));
					AssertThat(p_new.m_y, Equals(p_old.m_y+42.0));
				}
			}
		});

		it("scales",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.scale(-1.0, -2.0, true);
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(-ga.x(v)));
				AssertThat(grattr.y(v), Equals(-2.0 * ga.y(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(-p_old.m_x));
					AssertThat(p_new.m_y, Equals(-2.0 * p_old.m_y));
				}
			}
		});

		it("scales and then translates",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.scaleAndTranslate(-1.0, -42.0, 13, 37, true);
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(-ga.x(v) + 13));
				AssertThat(grattr.y(v), Equals(-42 * ga.y(v) + 37));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(-p_old.m_x + 13));
					AssertThat(p_new.m_y, Equals(-42 * p_old.m_y + 37));
				}
			}
		});

		it("flips vertical within its bounding box",[&](){
			DRect boundingBox = grattr.boundingBox();
			double height = boundingBox.height();
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.flipVertical();
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(ga.x(v)));
				AssertThat(grattr.y(v), Equals(height - ga.y(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(p_old.m_x));
					AssertThat(p_new.m_y, Equals(height - p_old.m_y));
				}
			}
		});

		it("flips vertical with a given box",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.flipVertical(DRect());
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(ga.x(v)));
				AssertThat(grattr.y(v), Equals(-ga.y(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(p_old.m_x));
					AssertThat(p_new.m_y, Equals(-p_old.m_y));
				}
			}
		});

		it("flips horizontal within its bounding box",[&](){
			DRect boundingBox = grattr.boundingBox();
			double width = boundingBox.width();
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.flipHorizontal();
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(width - ga.x(v)));
				AssertThat(grattr.y(v), Equals(ga.y(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(width - p_old.m_x));
					AssertThat(p_new.m_y, Equals(p_old.m_y));
				}
			}
		});

		it("flips horizontal with a given box",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.flipHorizontal(DRect());
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(-ga.x(v)));
				AssertThat(grattr.y(v), Equals(ga.y(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(-p_old.m_x));
					AssertThat(p_new.m_y, Equals(p_old.m_y));
				}
			}
		});

		it("rotates left",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.rotateLeft90();
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(ga.y(v)));
				AssertThat(grattr.y(v), Equals(-ga.x(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(p_old.m_y));
					AssertThat(p_new.m_y, Equals(-p_old.m_x));
				}
			}
		});

		it("rotates right",[&](){
			GraphAttributes ga = GraphAttributes(grattr);
			grattr.rotateRight90();
			for(node v : graph.nodes){
				AssertThat(grattr.x(v), Equals(-ga.y(v)));
				AssertThat(grattr.y(v), Equals(ga.x(v)));
			}
			for(edge e : graph.edges){
				DPolyline &bendpoints = ga.bends(e);
				for(DPoint &p_new : grattr.bends(e)){
					DPoint p_old = bendpoints.popFrontRet();
					AssertThat(p_new.m_x, Equals(-p_old.m_y));
					AssertThat(p_new.m_y, Equals(p_old.m_x));
				}
			}
		});
	});

	it("knows its bounding box",[&](){
		graph = Graph();
		randomGraph(graph, 100, 1000);
		grattr = GraphAttributes(graph, 0);
		grattr.initAttributes(Type::nodeGraphics | Type::edgeGraphics);
		for(node v : graph.nodes){
			grattr.x(v) = ogdf::randomNumber(-1000, 1000);
			grattr.y(v) = ogdf::randomNumber(-1000, 1000);
		}
		grattr.addNodeCenter2Bends(1);
		grattr.translateToNonNeg();
		DRect boundBox = grattr.boundingBox();
		AssertThat(boundBox.p1().m_x, Equals(0) || IsGreaterThan(0));
		AssertThat(boundBox.p1().m_y, Equals(0) || IsGreaterThan(0));
		AssertThat(boundBox.p2().m_x, Equals(2020) || IsLessThan(2020));
		AssertThat(boundBox.p2().m_y, Equals(2020) || IsLessThan(2020));
		for(node v : graph.nodes){
			AssertThat(boundBox.contains(DPoint(grattr.x(v), grattr.y(v))), IsTrue());
		}
		for(edge e : graph.edges){
			for(DPoint &p : grattr.bends(e)){
				AssertThat(boundBox.contains(p), IsTrue());
			}
		}
	});

	describe("bends",[&](){
		before_each([&](){
			graph = Graph();
			completeGraph(graph, 3);
			grattr = GraphAttributes(graph, 0);
			grattr.initAttributes(Type::nodeGraphics | Type::edgeGraphics);
		});

		it("clears all bends",[&](){
			grattr.addNodeCenter2Bends(1);
			AssertThat(grattr.bends(graph.chooseEdge()).size(), IsGreaterThan(0));
			grattr.clearAllBends();
			for(edge e : graph.edges){
				AssertThat(grattr.bends(e).size(), Equals(0));
			}
		});

		it("knows its bends",[&](){
			for(edge e : graph.edges){
				AssertThat(grattr.bends(e).size(), Equals(0));
			}
			grattr.addNodeCenter2Bends(0);
			for(edge e : graph.edges){
				AssertThat(grattr.bends(e).size(), Equals(2));
			}
		});
	});
});
});
