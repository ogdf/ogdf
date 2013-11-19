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

#include <ogdf/fileformats/DotParser.h>
#include <ogdf/fileformats/Utils.h>


namespace ogdf {

namespace dot {


Ast::Graph::Graph(
	const bool &strict,
	const bool &directed,
	std::string *id,
	StmtList *statements)
: strict(strict), directed(directed), id(id), statements(statements)
{
}


Ast::Graph::~Graph()
{
	delete id;
	delete statements;
}


Ast::StmtList::StmtList(
	Stmt *head,
	StmtList *tail)
: head(head), tail(tail)
{
}


Ast::StmtList::~StmtList()
{
	delete head;
	delete tail;
}


Ast::Stmt::~Stmt()
{
}


Ast::NodeStmt::NodeStmt(
	NodeId *nodeId,
	AttrList *attrs)
: nodeId(nodeId), attrs(attrs)
{
}


Ast::NodeStmt::~NodeStmt()
{
	delete nodeId;
	delete attrs;
}


Ast::EdgeStmt::EdgeStmt(
	EdgeLhs *lhs,
	EdgeRhs *rhs,
	AttrList *attrs)
: lhs(lhs), rhs(rhs), attrs(attrs)
{
}


Ast::EdgeStmt::~EdgeStmt()
{
	delete lhs;
	delete rhs;
	delete attrs;
}


Ast::AsgnStmt::AsgnStmt(
	const std::string &lhs,
	const std::string &rhs)
: lhs(lhs), rhs(rhs)
{
}


Ast::AsgnStmt::~AsgnStmt()
{
}


Ast::AttrStmt::AttrStmt(
	const Type &type,
	AttrList *attrs)
: type(type), attrs(attrs)
{
}


Ast::AttrStmt::~AttrStmt()
{
	delete attrs;
}


Ast::Subgraph::Subgraph(
	std::string *id,
	StmtList *statements)
: id(id), statements(statements)
{
}


Ast::Subgraph::~Subgraph()
{
	delete id;
	delete statements;
}


Ast::EdgeLhs::~EdgeLhs()
{
}


Ast::EdgeRhs::EdgeRhs(
	EdgeLhs *head,
	EdgeRhs *tail)
: head(head), tail(tail)
{
}


Ast::EdgeRhs::~EdgeRhs()
{
	delete head;
	delete tail;
}


Ast::NodeId::NodeId(
	const std::string &id,
	Port *port)
: id(id), port(port)
{
}


Ast::NodeId::~NodeId()
{
	delete port;
}


Ast::Port::Port(
	std::string *id,
	CompassPt *compassPt)
: id(id), compassPt(compassPt)
{
}


Ast::Port::~Port()
{
	delete id;
	delete compassPt;
}


Ast::CompassPt::CompassPt(
	const Type &type)
: type(type)
{
}


Ast::CompassPt::~CompassPt()
{
}


Ast::AttrList::AttrList(
	AList *head,
	AttrList *tail)
: head(head), tail(tail)
{
}


Ast::AttrList::~AttrList()
{
	delete head;
	delete tail;
}


Ast::AList::AList(
	AsgnStmt *head,
	AList *tail)
: head(head), tail(tail)
{
}


Ast::AList::~AList()
{
	delete head;
	delete tail;
}


Ast::Ast(const Tokens &tokens)
: m_tokens(tokens), m_tend(tokens.end()), m_graph(NULL)
{
}


Ast::~Ast()
{
	delete m_graph;
}


bool Ast::build()
{
	Iterator it = m_tokens.begin();
	delete m_graph;
	return (m_graph = parseGraph(it, it)) != 0;
}


Ast::Graph *Ast::root() const
{
	return m_graph;
}


Ast::EdgeStmt *Ast::parseEdgeStmt(
	Iterator curr, Iterator &rest)
{
	EdgeLhs *lhs;
	if(!((lhs = parseNodeId(curr, curr)) ||
	     (lhs = parseSubgraph(curr, curr)))) {
		return NULL;
	}

	EdgeRhs *rhs = parseEdgeRhs(curr, curr);
	if(!rhs) {
		delete lhs;
		return NULL;
	}

	AttrList *attrs = parseAttrList(curr, curr);

	rest = curr;
	return new EdgeStmt(lhs, rhs, attrs);
}


Ast::EdgeRhs *Ast::parseEdgeRhs(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || (curr->type != Token::edgeOpDirected &&
		                  curr->type != Token::edgeOpUndirected)) {
		return NULL;
	}
	curr++;

	EdgeLhs *head;
	if(!((head = parseSubgraph(curr, curr)) ||
	     (head = parseNodeId(curr, curr)))) {
		return NULL;
	}

	EdgeRhs *tail = parseEdgeRhs(curr, curr);

	rest = curr;
	return new EdgeRhs(head, tail);
}


Ast::NodeStmt *Ast::parseNodeStmt(
	Iterator curr, Iterator &rest)
{
	NodeId *nodeId = parseNodeId(curr, curr);
	if(!nodeId) {
		return NULL;
	}

	AttrList *attrs = parseAttrList(curr, curr);

	rest = curr;
	return new NodeStmt(nodeId, attrs);
}


Ast::NodeId *Ast::parseNodeId(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || curr->type != Token::identifier) {
		return NULL;
	}
	std::string id = *(curr->value);
	curr++;

	Port *port = parsePort(curr, curr);

	rest = curr;
	return new NodeId(id, port);
}


Ast::CompassPt *Ast::parseCompassPt(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || curr->type != Token::identifier) {
		return NULL;
	}
	const std::string &str = *(curr->value);
	curr++;
	if(str == "n") {
		curr = rest;
		return new CompassPt(CompassPt::n);
	}
	if(str == "ne") {
		rest = curr;
		return new CompassPt(CompassPt::ne);
	}
	if(str == "e") {
		rest = curr;
		return new CompassPt(CompassPt::e);
	}
	if(str == "se") {
		rest = curr;
		return new CompassPt(CompassPt::se);
	}
	if(str == "s") {
		rest = curr;
		return new CompassPt(CompassPt::s);
	}
	if(str == "sw") {
		rest = curr;
		return new CompassPt(CompassPt::sw);
	}
	if(str == "w") {
		rest = curr;
		return new CompassPt(CompassPt::w);
	}
	if(str == "nw") {
		rest = curr;
		return new CompassPt(CompassPt::nw);
	}
	if(str == "c") {
		rest = curr;
		return new CompassPt(CompassPt::c);
	}
	if(str == "_") {
		rest = curr;
		return new CompassPt(CompassPt::wildcard);
	}
	return NULL;
}


Ast::Port *Ast::parsePort(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || curr->type != Token::colon) {
		return NULL;
	}
	curr++;

	CompassPt *compass = parseCompassPt(curr, curr);
	if(compass) {
		rest = curr;
		return new Port(NULL, compass);
	}

	std::string *id = curr->value;
	curr++;

	if(curr != m_tend && curr->type == Token::colon) {
		curr++;

		compass = parseCompassPt(curr, curr);
		if(compass) {
			rest = curr;
			return new Port(id, compass);
		}

		// Compass parsing not succeeded, "put back" the colon.
		curr--;
	}

	rest = curr;
	return new Port(id, NULL);
}


Ast::AttrStmt *Ast::parseAttrStmt(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend) {
		return NULL;
	}

	AttrStmt::Type type;
	switch(curr->type) {
	case Token::graph:
		type = AttrStmt::graph;
		break;
	case Token::node:
		type = AttrStmt::node;
		break;
	case Token::edge:
		type = AttrStmt::edge;
		break;
	default:
		return NULL;
	}
	curr++;

	AttrList *attrs = parseAttrList(curr, curr);
	if(!attrs) {
		return NULL;
	}

	rest = curr;
	return new AttrStmt(type, attrs);
}


Ast::AsgnStmt *Ast::parseAsgnStmt(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || curr->type != Token::identifier) {
		return NULL;
	}
	std::string lhs = *(curr->value);
	curr++;

	if(curr == m_tend || curr->type != Token::assignment) {
		return NULL;
	}
	curr++;

	if(curr == m_tend || curr->type != Token::identifier) {
		return NULL;
	}
	std::string rhs = *(curr->value);
	curr++;

	rest = curr;
	return new AsgnStmt(lhs, rhs);
}


Ast::Subgraph *Ast::parseSubgraph(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend) {
		return NULL;
	}

	// Optional "subgraph" keyword and optional identifier.
	std::string *id = NULL;
	if(curr->type == Token::subgraph) {
		curr++;
		if(curr == m_tend) {
			return NULL;
		}
		if(curr->type == Token::identifier) {
			id = new std::string(*(curr->value));
			curr++;
		}
	}

	if(curr == m_tend || curr->type != Token::leftBrace) {
		delete id;
		return NULL;
	}
	curr++;

	StmtList *stmts = parseStmtList(curr, curr);

	if(curr == m_tend || curr->type != Token::rightBrace) {
		delete id;
		delete stmts;
		return NULL;
	}
	curr++;

	rest = curr;
	return new Subgraph(id, stmts);;
}


Ast::Stmt *Ast::parseStmt(
	Iterator curr, Iterator &rest)
{
	Stmt *stmt;
	if((stmt = parseEdgeStmt(curr, curr)) ||
	   (stmt = parseNodeStmt(curr, curr)) ||
	   (stmt = parseAttrStmt(curr, curr)) ||
	   (stmt = parseAsgnStmt(curr, curr)) ||
	   (stmt = parseSubgraph(curr, curr))) {
		rest = curr;
		return stmt;
	}

	return NULL;
}


Ast::StmtList *Ast::parseStmtList(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend) {
		return NULL;
	}

	Stmt *head = parseStmt(curr, curr);
	if(!head) {
		return NULL;
	}

	// Optional semicolon.
	if(curr != m_tend && curr->type == Token::semicolon) {
		curr++;
	}

	StmtList *tail = parseStmtList(curr, curr);
	rest = curr;
	return new StmtList(head, tail);
}


Ast::Graph *Ast::parseGraph(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend) {
		return NULL;
	}

	bool strict = false;
	bool directed = false;
	std::string *id = NULL;

	if(curr->type == dot::Token::strict) {
		strict = true;
		curr++;
	}

	if(curr == m_tend) {
		return NULL;
	}

	switch(curr->type) {
	case Token::graph:
		directed = false;
		break;
	case Token::digraph:
		directed = true;
		break;
	default:
		std::cerr << "ERROR: Unexpected token \""
		          << Token::toString(curr->type)
		          << "\" at "
		          << curr->row << ", " << curr->column << ".\n";
		return NULL;
	}
	curr++;

	if(curr == m_tend) {
		return NULL;
	}

	if(curr->type == Token::identifier) {
		id = new std::string(*(curr->value));
		curr++;
	}

	if(curr == m_tend || curr->type != Token::leftBrace) {
		// cerr << "ERROR: Expected \""
		//      << Token::toString(Token::leftBrace)
		//      << ", found \""
		//      << Token::toString(a->type)
		//      << "\" at "
		//      << a->row << ", " << a->column << ".\n";
		delete id;
		return NULL;
	}
	curr++;

	StmtList *statements = parseStmtList(curr, curr);

	if(curr == m_tend || curr->type != Token::rightBrace) {
		std::cerr << "ERROR: Expected \""
		          << Token::toString(Token::rightBrace)
		          << ", found \""
		          << Token::toString(curr->type)
		          << "\" at "
		          << curr->row << ", " << curr->column << ".\n";
		delete id;
		delete statements;
		return NULL;
	}
	curr++;

	rest = curr;
	return new Graph(strict, directed, id, statements);
}


Ast::AttrList *Ast::parseAttrList(
	Iterator curr, Iterator &rest)
{
	if(curr == m_tend || curr->type != Token::leftBracket) {
		return NULL;
	}
	curr++;

	AList *head = parseAList(curr, curr);

	if(curr == m_tend || curr->type != Token::rightBracket) {
		delete head;
		return NULL;
	}
	curr++;

	AttrList *tail = parseAttrList(curr, curr);

	rest = curr;
	return new AttrList(head, tail);
}


Ast::AList *Ast::parseAList(
	Iterator curr, Iterator &rest)
{
	AsgnStmt *head = parseAsgnStmt(curr, curr);
	if(!head) {
		return NULL;
	}

	// Optional comma.
	if(curr != m_tend && curr->type == Token::comma) {
		curr++;
	}

	AList *tail = parseAList(curr, curr);

	rest = curr;
	return new AList(head, tail);
}


static bool readBends(
	const std::string &str,
	DPolyline &polyline)
{
	// First, just trim every unnecessary charater - we don't treat DOT's
	// spline as spline but just set of bending points. One can always
	// implement B-splines and then generate bending points.
	std::string fixed(str);
	for(size_t i = 0; i < fixed.size(); i++) {
		if(fixed[i] == ',' || fixed[i] == ';' ||
		   fixed[i] == 'e' || fixed[i] == 'p')
		{
			fixed[i] = ' ';
		}
	}

	std::istringstream is(fixed);

	double x, y;
	polyline.clear();
	while(is >> x && is >> y) {
		polyline.pushBack(DPoint(x, y));
	}

	return true;
}


static bool readAttribute(
	GraphAttributes &GA, const node &v,
	const Ast::AsgnStmt &stmt)
{
	const long flags = GA.attributes();

	std::istringstream ss(stmt.rhs);
	switch(toAttribute(stmt.lhs)) {
	case a_label:
		if(flags & GraphAttributes::nodeLabel) {
			GA.label(v) = stmt.rhs;
		}
		break;
	case a_template:
		if(flags & GraphAttributes::nodeTemplate) {
			GA.templateNode(v) = stmt.rhs;
		}
		break;
	case a_width:
		if(flags & GraphAttributes::nodeGraphics) {
			// sscanf(stmt.rhs.c_str(), "%lf", &GA.width(v));
			ss >> GA.width(v);
		}
		break;
	case a_height:
		if(flags & GraphAttributes::nodeGraphics) {
			// sscanf(stmt.rhs.c_str(), "%lf", &GA.height(v));
			ss >> GA.height(v);
		}
		break;
	case a_shape:
		if(flags & GraphAttributes::nodeGraphics) {
			GA.shape(v) = toShape(stmt.rhs);
		}
		break;
	case a_position:
		if(flags & GraphAttributes::nodeGraphics) {
			// sscanf(stmt.rhs.c_str(), "%lf,%lf", &GA.x(v), &GA.y(v));
			ss >> GA.x(v) >> ',' >> GA.y(v);
		}
		break;
	case a_stroke:
		if(flags & GraphAttributes::nodeStyle) {
			GA.strokeColor(v) = stmt.rhs;
			// TODO: color literals.
		}
		break;
	case a_fill:
		if(flags & GraphAttributes::nodeStyle) {
			GA.fillColor(v) = stmt.rhs;
			// TODO: color literals.
		}
		break;
	default:
		std::cerr << "WARNING: Attribute \"" << stmt.lhs
		          << "\" is  not supported by node or incorrect. Ignoring.\n";
	}

	return true;
}


static bool readAttribute(
	GraphAttributes &GA, const edge &e,
	const Ast::AsgnStmt &stmt)
{
	const long flags = GA.attributes();

	std::istringstream ss(stmt.rhs);
	switch(toAttribute(stmt.lhs)) {
	case a_label:
		if(flags & GraphAttributes::edgeLabel) {
			GA.label(e) = stmt.rhs;;
		}
		break;
	case a_weight:
		if(flags & GraphAttributes::edgeDoubleWeight) {
			// sscanf(stmt.rhs.c_str(), "%lf", &GA.doubleWeight(e));
			ss >> GA.doubleWeight(e);
		} else if(flags & GraphAttributes::edgeIntWeight) {
			// sscanf(stmt.rhs.c_str(), "%d", &GA.intWeight(e));
			ss >> GA.intWeight(e);
		}
		break;
	case a_position:
		if(flags & GraphAttributes::edgeGraphics) {
			readBends(stmt.rhs, GA.bends(e));
		}
		break;
	case a_stroke:
		if(flags & GraphAttributes::edgeStyle) {
			GA.strokeColor(e) = stmt.rhs;
			// TODO: color literals.
		}
		break;
	case a_arrow:
		if(flags & GraphAttributes::edgeArrow) {
			GA.arrowType(e) = toArrow(stmt.rhs);
		}
	default:
		std::cerr << "WARNING: Attribute \"" << stmt.lhs
		          << "\" is not supported by edge or incorrect. Ignoring.\n";
	}

	return true;
}


static bool readAttribute(
	ClusterGraphAttributes &CA, const cluster &c,
	const Ast::AsgnStmt &stmt)
{
	switch(toAttribute(stmt.lhs)) {
	case a_label:
		CA.label(c) = stmt.rhs;
		break;
	case a_template:
		CA.templateCluster(c) = stmt.rhs;
		break;
	case a_fill:
		CA.fillColor(c) = stmt.rhs;
		break;
	case a_stroke:
		CA.strokeColor(c) = stmt.rhs;
		break;
	default:
		std::cerr << "WARNING: Attribute \"" << stmt.lhs
	              << "\" is not supported by cluster or incorrect. Ignoring.\n";
	}
	return true;
}


template <typename G, typename T>
static inline bool readAttributes(
	G &GA, T elem,
	Ast::AttrList *attrs)
{
	for(; attrs; attrs = attrs->tail) {
		for(Ast::AList *alist = attrs->head; alist; alist = alist->tail) {
			if(!readAttribute(GA, elem, *(alist->head))) {
				return false;
			}
		}
	}

	return true;
}


template <typename G, typename T>
static inline bool readAttributes(
	G &GA, T elem,
	const std::vector<Ast::AttrList *> &defaults)
{
	for(std::vector<Ast::AttrList *>::const_iterator it = defaults.begin();
	    it != defaults.end();
	    it++)
	{
		if(!readAttributes(GA, elem, *it)) {
			return false;
		}
	}

	return true;
}


static inline bool readStatements(
	Parser &P,
	Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data,
	Ast::StmtList *stmts)
{
	for(; stmts; stmts = stmts->tail) {
		if(!stmts->head->read(P, G, GA, C, CA, data)) {
			return false;
		}
	}

	return true;
}


bool Ast::Graph::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA)
{
	if(GA) {
		GA->setDirected(directed);
	}

	std::set<node> subgraphNodes;
	std::vector<AttrList *> nodeDefaults, edgeDefaults;
	return readStatements(
		P, G, GA, C, CA,
		SubgraphData(
			// Root cluster.
			C ? C->rootCluster() : NULL,
			nodeDefaults, edgeDefaults, subgraphNodes),
		statements);
}


bool Ast::NodeStmt::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	const node v = P.requestNode(G, GA, C, data, nodeId->id);
	data.nodes.insert(v);
	return GA ? readAttributes(*GA, v, attrs) : true;
}


static inline bool cross(
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const std::vector<Ast::AttrList *> &defaults, Ast::AttrList *attrs,
	const std::set<ogdf::node> &lnodes, const std::set<ogdf::node> &rnodes)
{
	for(std::set<node>::const_iterator it = lnodes.begin();
	    it != lnodes.end();
	    it++)
	{
		for(std::set<node>::const_iterator jt = rnodes.begin();
		    jt != rnodes.end();
		    jt++)
		{
			const edge e = G.newEdge(*it, *jt);
			if(GA && !(readAttributes(*GA, e, defaults) &&
			           readAttributes(*GA, e, attrs))) {
				return false;
			}
		}
	}

	return true;
}

bool Ast::EdgeStmt::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	Ast::EdgeLhs *lhs = this->lhs;

	std::set<node> lnodes;
	lhs->read(P, G, GA, C, CA, data.withNodes(lnodes));

	for(Ast::EdgeRhs *rhs = this->rhs; rhs; rhs = rhs->tail) {
		std::set<node> rnodes;
		rhs->head->read(P, G, GA, C, CA, data.withNodes(rnodes));

		if(!cross(G, GA, C, CA, data.edgeDefaults, attrs, lnodes, rnodes)) {
			return false;
		}

		// Append left side nodes to the result and make right node left ones.
		data.nodes.insert(lnodes.begin(), lnodes.end());
		std::swap(lnodes, rnodes);
		lhs = rhs->head;
	}

	return true;
}


bool Ast::AsgnStmt::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	return CA ? readAttribute(*CA, data.rootCluster, *this) : true;
}


bool Ast::AttrStmt::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	switch(type) {
	case graph:
		return CA ? readAttributes(*CA, data.rootCluster, attrs) : true;
	case node:
		data.nodeDefaults.push_back(attrs);
		return true;
	case edge:
		data.edgeDefaults.push_back(attrs);
		return true;
	default:
		return false;
	}
}


bool Ast::Subgraph::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	// We pass a copy of defaults to subgraph because those parameters should
	// be local for given subgraph.
	std::vector<AttrList *> nodeDefaults(data.nodeDefaults);
	std::vector<AttrList *> edgeDefaults(data.edgeDefaults);
	SubgraphData newData = data.withDefaults(nodeDefaults, edgeDefaults);

	// Create new cluster if subgraph identifier is given and it starts with
	// pattern "cluster". Otherwise, subgraph is not considered as a cluster
	// (as stated in DOT manual).
	const std::string patt = "cluster";
	if(C && id && id->compare(0, patt.length(), patt) == 0) {
		return readStatements(
			P, G, GA, C, CA,
			newData.withCluster(C->newCluster(newData.rootCluster)),
			statements);
	}

	return readStatements(P, G, GA, C, CA, newData, statements);
}


bool Ast::NodeId::read(
	Parser &P,
	ogdf::Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA,
	const SubgraphData &data)
{
	data.nodes.insert(P.requestNode(G, GA, C, data, id));
	return true;
}


Parser::Parser(std::istream &in) : m_in(in)
{
	m_nodeId = HashArray<std::string, node>(NULL);
}


node Parser::requestNode(
	Graph &G, GraphAttributes *GA, ClusterGraph *C,
	const SubgraphData &data,
	const std::string &id)
{
	node v;
	// Ugly and slow, fix it somehow in the future.
	if(!m_nodeId[id]) {
		v = m_nodeId[id] = G.newNode();
		if(C) {
			C->reassignNode(v, data.rootCluster);
		}

		// We also set default attributes when a node is requested for the
		// first time. This may sound strange but Graphviz's DOT tool behaves
		// exactly the same so I guess this is a right place to use these.
		if(GA) {
			if(GA->attributes() & GraphAttributes::nodeLabel) {
				GA->label(v) = id;
			}
			readAttributes(*GA, v, data.nodeDefaults);
		}
	} else {
		v = m_nodeId[id];
	}

	// So, the question is: where to put a node if it can be declared with
	// edge statement and edge statements with common node may appear in
	// various clusters? Accoring to my own tests (using Graphviz's DOT tool)
	// it is simple: put it in the deepest possible cluster in which it shows
	// up. And this is achieved by the line below - if a node is requested on
	// a level that is deeper than currently assigned one, then we reassign
	// cluster.
	if(C && data.rootCluster->depth() < C->clusterOf(v)->depth()) {
		C->reassignNode(v, data.rootCluster);
	}

	return v;
}


bool Parser::readGraph(
	Graph &G, GraphAttributes *GA,
	ClusterGraph *C, ClusterGraphAttributes *CA)
{
	m_nodeId.clear();
	G.clear();
	if(C) {
		C->semiClear();
	}

	Lexer lexer(m_in);
	if(!lexer.tokenize()) {
		return false;
	}

	Ast ast(lexer.tokens());
	if(!ast.build()) {
		return false;
	}

	return ast.root()->read(*this, G, GA, C, CA);
}


bool Parser::read(Graph &G)
{
	return readGraph(G, NULL, NULL, NULL);
}


bool Parser::read(Graph &G, GraphAttributes &GA)
{
	return readGraph(G, &GA, NULL, NULL);
}


bool Parser::read(Graph &G, ClusterGraph &C)
{
	return readGraph(G, NULL, &C, NULL);
}


bool Parser::read(Graph &G, ClusterGraph &C, ClusterGraphAttributes &CA)
{
	return readGraph(G, &CA, &C, &CA);
}


SubgraphData::SubgraphData(
	cluster rootCluster,
	std::vector<Ast::AttrList *> &nodeDefaults,
	std::vector<Ast::AttrList *> &edgeDefaults,
	std::set<node> &nodes)
:
	rootCluster(rootCluster),
	nodeDefaults(nodeDefaults),
	edgeDefaults(edgeDefaults),
	nodes(nodes)
{
}


SubgraphData SubgraphData::withCluster(
	cluster newRootCluster) const
{
	return SubgraphData(newRootCluster, nodeDefaults, edgeDefaults, nodes);
}


SubgraphData SubgraphData::withDefaults(
	std::vector<Ast::AttrList *> &newNodeDefaults,
	std::vector<Ast::AttrList *> &newEdgeDefaults) const
{

	return SubgraphData(rootCluster, newNodeDefaults, newEdgeDefaults, nodes);
}


SubgraphData SubgraphData::withNodes(
	std::set<node> &newNodes) const
{
	return SubgraphData(rootCluster, nodeDefaults, edgeDefaults, newNodes);
}


} // end namespace dot

} // end namespace ogdf
