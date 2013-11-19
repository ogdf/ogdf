/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declarations for DOT Parser
 *
 * \author Łukasz Hanuszczak
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_DOT_PARSER_H
#define OGDF_DOT_PARSER_H

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/basic/HashArray.h>

#include <ogdf/fileformats/DOT.h>
#include <ogdf/fileformats/DotLexer.h>

#include <cstdio>
#include <vector>
#include <set>
#include <string>


namespace ogdf {

namespace dot {


class Ast;
class Parser;
struct SubgraphData;

//! DOT format abstract syntax tree class, based on official documentation.
/**
 * To provide easier modification and improve code readability the DOT format
 * parsing process is divided to three stages. Building an abstract syntax tree
 * is the second one. AST elements are also involved in reading to Graph,
 * GraphAttributes, ClusterGraph and ClusterGraphAttributes.
 *
 * This class provide method for building such tree from a token list. The AST
 * is then represented by tree of structures corresponding to DOT's language
 * abstract grammar. That grammar is nearly 100% like official one - minor
 * adjustments have been made just for my convenience. Neverthless, my version
 * should be equal to the official one in terms of represented language.
 *
<a href="http://www.graphviz.org/content/dot-language">
official DOT specification
</a>
 *
\verbatim
Graph = [ 'strict' ] ( 'graph' | 'digraph' ) [ id ] '{' [ StmtList ] '}'
Subgraph = [ 'subgraph' [ id ] ] '{' StmtList '}'

StmtList = Stmt [ ';' ] [ StmtList ]
Stmt = NodeStmt | EdgeStmt | AttrStmt | AsgnStmt | Subgraph

NodeStmt = NodeId [ AttrList ]
NodeId = id [ port ]

EdgeStmt = ( NodeId | Subgraph ) EdgeRhs [ AttrList ]
EdgeRhs = ( '--' | '->' ) ( NodeId | Subgraph) [ EdgeRhs ]

AttrStmt = ( 'graph' | 'node' | 'edge' ) AttrList

AsgnStmt = id '=' id

AttrList = '[' [ AList ] ']' [ AttrList ]
AList = AsgnStmt [ ',' ] [ AList ]

Port = ':' id [ ':' CompassPt ] | ':' CompassPt
CompassPt = ( 'n' | 'ne' | 'e' | 'se' | 's' | 'sw' | 'w' | 'nw' | 'c' | '_' )
\endverbatim
 *
 * The AST building process tries to mirror given grammar as much as possible
 * keeping the code clean and simple. Each grammar's nonterminal symbol has
 * corresponding \a parse function. These \a parse functions take two args:
 * \a current iterator and \a rest iterator. \a current iterator indicates
 * a position where the parsing should begin. On successfull parse, function
 * returns (pointer to) tree element and moves \a rest iterator to a place
 * where it has ended. On failure, function returns \c NULL pointer and does
 * nothing to \a rest iterator.
 *
 * Finally, non-list AST elements provide \a read methods. These functions
 * allow to, surprisingly, read Graph structure and/or associated attributes.
 * The entry point is Ast::Graph#read method which reads basic attributes
 * like graph strictness and triggers recursively \a read of its ancestors.
 *
 * \sa dot::Lexer
 * \sa dot::Parser
 */
class Ast {
public:
	struct Graph;
	struct StmtList;
	struct NodeStmt;
	struct EdgeStmt;
	struct AttrStmt;
	struct AsgnStmt;
	struct Subgraph;
	struct NodeId;
	struct EdgeRhs;
	struct AttrList;
	struct AList;
	struct Port;
	struct CompassPt;

	struct Stmt;
	struct EdgeLhs;

private:
	typedef std::vector<Token> Tokens;
	typedef Tokens::const_iterator Iterator;

	const Tokens m_tokens;
	const Iterator m_tend;

	Graph *m_graph;

	Graph *parseGraph(
		Iterator current, Iterator &rest);
	Subgraph *parseSubgraph(
		Iterator current, Iterator &rest);
	NodeStmt *parseNodeStmt(
		Iterator current, Iterator &rest);
	EdgeStmt *parseEdgeStmt(
		Iterator current, Iterator &rest);
	AttrStmt *parseAttrStmt(
		Iterator current, Iterator &rest);
	AsgnStmt *parseAsgnStmt(
		Iterator current, Iterator &rest);
	EdgeRhs *parseEdgeRhs(
		Iterator current, Iterator &rest);
	NodeId *parseNodeId(
		Iterator current, Iterator &rest);
	Stmt *parseStmt(
		Iterator current, Iterator &rest);
	StmtList *parseStmtList(
		Iterator current, Iterator &rest);
	AttrList *parseAttrList(
		Iterator current, Iterator &rest);
	AList *parseAList(
		Iterator current, Iterator &rest);
	Port *parsePort(
		Iterator current, Iterator &rest);
	CompassPt *parseCompassPt(
		Iterator current, Iterator &rest);

public:
	//! Initializes AST building but does not trigger the process itself.
	/**
	 * @param tokens DOT format token list to build the AST.
	 */
	Ast(const Tokens &tokens);
	~Ast();

	//! Builds the DOT format AST.
	/**
	 * @return True if success, false otherwise.
	 */
	bool build();

	//! Returns the root of the AST (NULL if none).
	Graph *root() const;

	struct Graph {
		const bool strict;
		const bool directed;
		std::string *id;

		StmtList *statements;

		Graph(
			const bool &strict,
			const bool &directed,
			std::string *id,
			StmtList *statements);
		~Graph();

		bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA);
	};

	struct StmtList {
		Stmt *head;
		StmtList *tail;

		StmtList(
			Stmt *head,
			StmtList *tail);
		~StmtList();
	};

	struct Stmt {
		virtual ~Stmt() = 0;

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data) = 0;
	};

	struct NodeStmt : public Stmt {
		NodeId *nodeId;
		AttrList *attrs;

		NodeStmt(
			NodeId *nodeId,
			AttrList *attrs);
		~NodeStmt();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct EdgeStmt : public Stmt {
		EdgeLhs *lhs;
		EdgeRhs *rhs;
		AttrList *attrs;

		EdgeStmt(
			EdgeLhs *lhs,
			EdgeRhs *rhs,
			AttrList *attrs);
		~EdgeStmt();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct AsgnStmt : public Stmt {
		const std::string lhs;
		const std::string rhs;

		AsgnStmt(
			const std::string &lhs,
			const std::string &rhs);
		~AsgnStmt();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct AttrStmt : public Stmt {
		enum Type { graph, edge, node };

		Type type;
		AttrList *attrs;

		AttrStmt(
			const Type &type,
			AttrList *attrs);
		~AttrStmt();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct EdgeLhs {
		virtual ~EdgeLhs() = 0;

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data) = 0;
	};

	struct Subgraph : public Stmt, EdgeLhs {
		std::string *id;
		StmtList *statements;

		Subgraph(
			std::string *id,
			StmtList *statements);
		~Subgraph();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct NodeId : public EdgeLhs {
		const std::string id;
		Port *port;

		NodeId(
			const std::string &id,
			Port *port);
		~NodeId();

		virtual bool read(
			Parser &P,
			ogdf::Graph &G, GraphAttributes *GA,
			ClusterGraph *C, ClusterGraphAttributes *CA,
			const SubgraphData &data);
	};

	struct CompassPt {
		enum Type { n, ne, e, se, s, sw, w, nw, c, wildcard };
		Type type;

		CompassPt(
			const Type &type);
		~CompassPt();
	};

	struct Port {
		std::string *id;
		CompassPt *compassPt;

		Port(
			std::string *id,
			CompassPt *compassPt);
		~Port();
	};

	struct EdgeRhs {
		EdgeLhs *head;
		EdgeRhs *tail;

		EdgeRhs(
			EdgeLhs *head,
			EdgeRhs *tail);
		~EdgeRhs();
	};

	struct AttrList {
		AList *head;
		AttrList *tail;

		AttrList(
			AList *head,
			AttrList *tail);
		~AttrList();
	};

	struct AList {
		AsgnStmt *head;
		AList *tail;

		AList(
			AsgnStmt *head,
			AList *tail);
		~AList();
	};
};


//! DOT format parser class.
/**
 * Provides methods for reading graphs in the DOT format. This class is
 *  actually just a wrapper/composer for Ast and Lexer classes.
 *
 * \sa dot::Lexer
 * \sa dot::Ast
 */
class Parser {
private:
	std::istream &m_in;

	// Maps node id to Graph node.
	HashArray<std::string, node> m_nodeId;

	bool readGraph(
		Graph &G, GraphAttributes *GA,
		ClusterGraph *C, ClusterGraphAttributes *CA);

public:
	//! Initializes parser class with given input (but does nothing to it).
	Parser(std::istream &in);

	bool read(Graph &G);
	bool read(Graph &G, GraphAttributes &GA);
	bool read(Graph &G, ClusterGraph &C);
	bool read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA);

	//! Perfoms a node \a query, returning node for given attribute.
	/**
	 * Returns a node with given id in a graph. If node is requested for the
	 * first time then Graph#newNode is called and node is initialized with
	 * default attributes and placed in proper cluster (through \a data).
	 * @param G Graph whom node is requested.
	 * @param GA GraphAttributes for given graph, ignored if \c NULL.
	 * @param C ClusterGraph for given graph, ignored if \c NULL.
	 * @param data Data about current subgraph.
	 * @param id Identifier of requested node.
	 * @return Requested node.
	 */
	node requestNode(
		Graph &G, GraphAttributes *GA, ClusterGraph *C,
		const SubgraphData &data,
		const std::string &id);
};


//! A helper structure containing information for recursive graph reading.
struct SubgraphData {
	cluster rootCluster;
	std::vector<Ast::AttrList *> &nodeDefaults;
	std::vector<Ast::AttrList *> &edgeDefaults;
	std::set<node> &nodes;

	//! Initializes structure with given data.
	/**
	 * @param rootCluster Root cluster of current subgraph.
	 * @param nodeDefaults Node default attributes.
	 * @param edgeDefaults Edge default attributes.
	 * @param nodes Nodes in current subgraph.
	 */
	SubgraphData(
		cluster rootCluster,
		std::vector<Ast::AttrList *> &nodeDefaults,
		std::vector<Ast::AttrList *> &edgeDefaults,
		std::set<node> &nodes);


	//! Returns almost the same structure, but with root cluster.
	SubgraphData withCluster(
		cluster newRootCluster) const;

	//! Returns almost the same structure, but with new defaults.
	SubgraphData withDefaults(
		std::vector<Ast::AttrList *> &newNodeDefaults,
		std::vector<Ast::AttrList *> &newEdgeDefaults) const;
	//! Returns almost the same structure, but with new node list.
	SubgraphData withNodes(
		std::set<node> &newNodes) const;
};


} // end namespace dot

} // end namespace ogdf


#endif
