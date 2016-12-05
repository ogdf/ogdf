/** \file
 * \brief Implementation of the class UmlToGraphConverter
 *
 * \author Dino Ahr
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


#include <ogdf/fileformats/UmlToGraphConverter.h>
#include <ogdf/fileformats/GraphIO.h>


namespace ogdf {


	//
	// C o n s t r u c t o r
	//
	UmlToGraphConverter::UmlToGraphConverter(istream &is)
	{
		// Create parser and get reference to hash table
		m_xmlParser = new XmlParser(is);

		// Fill hash table of the parser with predefined info indices
		initializePredefinedInfoIndices();

		// Create the parse tree
		if (m_xmlParser->createParseTree() == false) {
			GraphIO::logger.lout() << "Could not create XML parse tree!" << endl;
			return; // parse error
		}

		// Create the uml model graph
		m_modelGraph = new UmlModelGraph();
		if (!createModelGraph(*m_modelGraph)) {
			GraphIO::logger.lout() << "Could not create UML model graph." << endl;
			return;
		}

		// Create the uml diagram graphs
		if (!createDiagramGraphs()){
			GraphIO::logger.lout() << "Could not create UML diagram graphs." << endl;
			return;
		}

		// Create the diagram graph in UMLGraph format
		if (!createDiagramGraphsInUMLGraphFormat(m_diagramGraphsInUMLGraphFormat)) {
			GraphIO::logger.lout() << "Could not create diagram graph in UML graph format." << endl;
			return;
		}
	} // UmlToGraphConverter

	//
	// D e s t r u c t o r
	//
	UmlToGraphConverter::~UmlToGraphConverter()
	{
		// Delete diagram graphs in UMLGraph format
		for (UMLGraph *umlg : m_diagramGraphsInUMLGraphFormat) {
			const Graph & associatedGraph = (const Graph &)(*umlg);
			delete umlg;
			delete &associatedGraph;
		}
		m_diagramGraphsInUMLGraphFormat.clear();


		// Delete diagram graphs
		for (UmlDiagramGraph *dg : m_diagramGraphs) {
			delete dg;
		}
		m_diagramGraphs.clear();

		// Destroy model graph
		delete m_modelGraph;

		// Destroy parser
		delete m_xmlParser;
	}

	//
	// i n i t i a l i z e P r e d e f i n e d I n f o I n d i c e s
	//
	void UmlToGraphConverter::initializePredefinedInfoIndices()
	{
		m_xmlParser->addNewHashElement("XMI",                        xmi);
		m_xmlParser->addNewHashElement("XMI.content",			     xmiContent);
		m_xmlParser->addNewHashElement("xmi.id",                     xmiId);
		m_xmlParser->addNewHashElement("UML:Model",                  umlModel);
		m_xmlParser->addNewHashElement("UML:Namespace.ownedElement", umlNamespaceOwnedElement);
		m_xmlParser->addNewHashElement("UML:Class",					 umlClass);
		m_xmlParser->addNewHashElement("name",						 name);
		m_xmlParser->addNewHashElement("UML:Generalization",		 umlGeneralization);
		m_xmlParser->addNewHashElement("child",						 child);
		m_xmlParser->addNewHashElement("parent",					 parent);
		m_xmlParser->addNewHashElement("UML:Association",			 umlAssociation);
		m_xmlParser->addNewHashElement("UML:Association.connection", umlAssociationConnection);
		m_xmlParser->addNewHashElement("UML:AssociationEnd",		 umlAssociationEnd);
		m_xmlParser->addNewHashElement("type",		                 type);
		m_xmlParser->addNewHashElement("UML:Diagram",		         umlDiagram);
		m_xmlParser->addNewHashElement("UML:Diagram.element",		 rootUmlDiagramElement);
		m_xmlParser->addNewHashElement("UML:DiagramElement",		 umlDiagramElement);
		m_xmlParser->addNewHashElement("geometry",					 geometry);
		m_xmlParser->addNewHashElement("subject",					 subject);
		m_xmlParser->addNewHashElement("UML:Package",				 umlPackage);
		m_xmlParser->addNewHashElement("UML:Interface",				 umlInterface);
		m_xmlParser->addNewHashElement("UML:Dependency",			 umlDependency);
		m_xmlParser->addNewHashElement("client",			         client);
		m_xmlParser->addNewHashElement("supplier",			         supplier);
		m_xmlParser->addNewHashElement("diagramType",			     diagramType);
		m_xmlParser->addNewHashElement("ClassDiagram",			     classDiagram);
		m_xmlParser->addNewHashElement("ModuleDiagram",			     moduleDiagram);

	} // initializePredefinedInfoIndices


	//
	// p r i n t I d T o N o d e M a p p i n g T a b l e
	//
	void UmlToGraphConverter::printIdToNodeMappingTable(ofstream &os)
	{
		// Header
		os << "\n--- Content of Hash table: m_m_idToNode ---\n" << endl;

		// Get iterator
		HashConstIterator<int, NodeElement*> it;

		// Traverse table
		for( it = m_idToNode.begin(); it.valid(); ++it){
			os << "\"" << it.key() << "\" has index "
				<< m_modelGraph->getNodeLabel(it.info()) << endl;
		}

	} // printIdToNodeMappingTable

	//
	// p r i n t D i a g r a m s I n U M L G r a p h F o r m a t
	//
	void UmlToGraphConverter::printDiagramsInUMLGraphFormat(ofstream &os)
	{
		// Traverse diagrams
		for (UMLGraph *diagram : m_diagramGraphsInUMLGraphFormat)
		{
			// Get underlying graphs
			const Graph &G = (const Graph &)*diagram;
			const GraphAttributes &AG = *diagram;

			// Nodes
			os << "Classes:" << endl;
			for(node v : G.nodes)
			{
				os << "\t" << AG.label(v);

				os << " with geometry ("
					 << AG.x(v) << ", "
					 << AG.y(v) << ", "
					 << AG.width(v) << ", "
					 << AG.height(v) << ")";

				os << endl;
			}

			// Edges
			os << "Relations:" << endl;
			for(edge e : G.edges)
			{
				os << "\t";

				if (AG.type(e) == Graph::association)
					os << "Association between ";
				if (AG.type(e) == Graph::generalization)
					os << "Generalization between ";

				os << AG.label(e->source()) << " and "
					 << AG.label(e->target()) << endl;
			}

			os << "---------------------------------------------------------------\n\n" << endl;

		} // Traverse diagrams

	} // printDiagramsInUMLGraphFormat


	//
	// c r e a t e M o d e l G r a p h
	//
	bool UmlToGraphConverter::createModelGraph(UmlModelGraph &modelGraph){

		// Message
		//cout << "Creating model graph..." << endl;

		// Check root element (must be <XMI>)
		if (m_xmlParser->getRootTag().m_pTagName->info() != xmi) {
			GraphIO::logger.lout() << "Root tag is not <XMI>" << endl;
			return false;
		}

		// Find first <UML:Namespace.ownedElement>; this is the father tag
		Array<int> path(3);
		path[0] = xmiContent;
		path[1] = umlModel;
		path[2] = umlNamespaceOwnedElement;
		const XmlTagObject *fatherTag;
		string rootPackageName("");
		if (!m_xmlParser->traversePath(m_xmlParser->getRootTag(), path, fatherTag)) {
			GraphIO::logger.lout() << "Path xmiContent, umlModel, umlNamespaceOwnedElement not found!" << endl;
			return false;
		}

		// Traverse packages and insert classifier nodes
		if (!traversePackagesAndInsertClassifierNodes(
			*fatherTag,
			rootPackageName,
			modelGraph))
		{
			return false;
		}

		// Note that first alle nodes have to be inserted into the model graph
		// and after that the edges should be inserted. The reason is that it
		// is possible that edges are specified prior to that one or both nodes
		// have been created.


		// Traverse packages and insert association edges
		if (!traversePackagesAndInsertAssociationEdges(*fatherTag, modelGraph))
		{
			return false;
		}

		// Traverse packages and insert generalization edges
		if (!traversePackagesAndInsertGeneralizationEdges(*fatherTag, modelGraph))
		{
			return false;
		}

		// Insert dependency edges
		if (!insertDependencyEdges(*fatherTag, modelGraph))
		{
			return false;
		}

		return true;

	} // createModelGraph


	//
	// t r a v e r s e P a c k a g e s A n d I n s e r t C l a s s i f i e r N o d e s
	//
	bool UmlToGraphConverter::traversePackagesAndInsertClassifierNodes(
		const XmlTagObject &currentRootTag,
		const string &currentPackageName,
		UmlModelGraph &modelGraph)
	{
		// We proceed in a DFS manner. As long as we are inside a package
		// and there is a subpackage inside, we dive into that subpackage
		// by calling this function recursively (with a new rootTag). Along
		// this we also construct the appropriate package name.
		//
		// If we arrive at a level where either all subpackages have been
		// already traversed or no subpackage is inside we proceed to find
		// the classifiers contained at the current level.
		//
		// At the moment we consider classed and interfaces.
		//
		// TODO: In Java it is possible that there classes contained in other classe,
		//       This is currently not detected.

		// Identify contained packages (<UML:Package>)
		const XmlTagObject *packageSon = nullptr;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlPackage, packageSon);
		while(packageSon != nullptr){

			// Create new name for the subpackage
			const XmlAttributeObject *nameAttribute;
			m_xmlParser->findXmlAttributeObject(*packageSon, name, nameAttribute);
			OGDF_ASSERT(nameAttribute != nullptr);
			string subPackageName = currentPackageName;
			if (currentPackageName.length() != 0){
				subPackageName += "::";
			}
			subPackageName += nameAttribute->m_pAttributeValue->key();

			// Find son umlNamespaceOwnedElement which indicates a nested package
			// if nonexistent then continue
			const XmlTagObject *newRootTag;
			if(m_xmlParser->findSonXmlTagObject(*packageSon, umlNamespaceOwnedElement, newRootTag)){

				// Call this function recursively
				if (!traversePackagesAndInsertClassifierNodes(*newRootTag, subPackageName, modelGraph))
				{
					// Something went wrong
					return false;
				}

			}

			// Next package (will be put into packageSon)
			m_xmlParser->findBrotherXmlTagObject(*packageSon, umlPackage, packageSon);

		} // while

		// Identify contained classes (<UML:Class>)
		if (!insertSpecificClassifierNodes(currentRootTag, currentPackageName, umlClass, modelGraph))
		{
			// Something went wrong
			return false;
		}

		// Identify contained interfaces (<UML:Interface>)
		if (!insertSpecificClassifierNodes(currentRootTag, currentPackageName, umlInterface, modelGraph))
		{
			// Something went wrong
			return false;
		}

		return true;

	} // traversePackagesAndInsertClassifierNodes


	//
	// i n s e r t S p e c i f i c C l a s s i f i e r N o d e s
	//
	bool UmlToGraphConverter::insertSpecificClassifierNodes(const XmlTagObject &currentRootTag,
																const string currentPackageName,
																int desiredClassifier,
																UmlModelGraph &modelGraph)
	{
		const XmlTagObject *classifierSon;
		m_xmlParser->findSonXmlTagObject(currentRootTag, desiredClassifier, classifierSon);
		while (classifierSon != nullptr){

			// Use the infoIndex of value of attribute xmi.id as reference to the node
			// it is unique for each classifier and is used to reference it in the
			// relation specifications
			const XmlAttributeObject *xmiIdAttr;

			// Did not find attribute xmi.id of classifier
			if (!m_xmlParser->findXmlAttributeObject(*classifierSon, xmiId, xmiIdAttr)) {
				GraphIO::logger.lout() << "Did not find attribute xmi.id of classifier." << endl;
				return false;
			}

			// We get an unique node id by the value of attribute xmi.id
			int nodeId = xmiIdAttr->m_pAttributeValue->info();

			// Find out name of the classifier
			const XmlAttributeObject *nameAttr;

			// Did not find name attribute
			if (!m_xmlParser->findXmlAttributeObject(*classifierSon, name, nameAttr)) {
				GraphIO::logger.lout() << "Did not find name attribute of classifier." << endl;
				return false;
			}

			// Name of the classifier is contained in the tag value
			HashedString *nodeName = nameAttr->m_pAttributeValue;

			// Create classifier name by prefixing it with the package name
			string nodeNameString = currentPackageName;
			if (currentPackageName.length() != 0){
				nodeNameString += "::";
			}
			nodeNameString += nodeName->key();

			// Check if node already exists
			if (m_idToNode.lookup(nodeId) != nullptr) {
				GraphIO::logger.lout() << "Node already exists." << endl;
				return false;
			}

			// Create a node for the graph
			NodeElement *node = modelGraph.newNode();
			modelGraph.label(node) = nodeNameString;
			modelGraph.type(node) = Graph::vertex;

			// Put node into hash table
			m_idToNode.fastInsert(nodeId, node);

			// Proceed with next class (will be put into classifierSon)
			m_xmlParser->findBrotherXmlTagObject(*classifierSon, desiredClassifier, classifierSon);

		} // while (classifierSon != 0)

		return true;

	} // insertSpecificClassifierNodes


	//
	// t r a v e r s e P a c k a g e s A n d I n s e r t A s s o c i a t i o n E d g e s
	//
	bool UmlToGraphConverter::traversePackagesAndInsertAssociationEdges(const XmlTagObject &currentRootTag,
																			UmlModelGraph &modelGraph)
	{
		// The traversion of the packages is identical with this of
		// traversePackagesAndInsertClassifierNodes

		// Identify contained packages (<UML:Package>)
		const XmlTagObject *packageSon;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlPackage, packageSon);
		while (packageSon != nullptr){

			// Find son umlNamespaceOwnedElement
			// if nonexistent then continue
			const XmlTagObject *newRootTag;

			if (m_xmlParser->findSonXmlTagObject(*packageSon, umlNamespaceOwnedElement, newRootTag))
			{
				// Call this function recursively
				if (!traversePackagesAndInsertAssociationEdges(*newRootTag, modelGraph))
				{
					return false;
				}

			}

			// Next package
			m_xmlParser->findBrotherXmlTagObject(*packageSon, umlPackage, packageSon);

		} // while

		// Find all associations (<UML:Association>)
		const XmlTagObject *associationSon;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlAssociation, associationSon);
		while (associationSon != nullptr){

			// Find out the reference number of this edge
			const XmlAttributeObject *edgeIdAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*associationSon, xmiId, edgeIdAttr);
			int edgeId = edgeIdAttr->m_pAttributeValue->info();

			// Go to <UML:Association.connection>
			const XmlTagObject *connection;
			m_xmlParser->findSonXmlTagObject(*associationSon, umlAssociationConnection, connection);

			// We assume binary associations

			// Investigate association ends
			const XmlTagObject *end1;
			m_xmlParser->findSonXmlTagObject(*connection, umlAssociationEnd, end1);

			// Something wrong
			if (!end1) {
				GraphIO::logger.lout(Logger::LL_MINOR) << "Current association tag does not contain both end tags!" << endl;
				// Next association
				m_xmlParser->findBrotherXmlTagObject(*associationSon, umlAssociation, associationSon);
				continue;
			}

			const XmlTagObject *end2;
			m_xmlParser->findBrotherXmlTagObject(*end1, umlAssociationEnd, end2);

			// Something wrong
			if (!end2) {
				GraphIO::logger.lout(Logger::LL_MINOR) << "Current association tag does not contain both end tags!" << endl;
				// Next association
				m_xmlParser->findBrotherXmlTagObject(*associationSon, umlAssociation, associationSon);
				continue;
			}

			// Use the infoIndex of value of attribute type to find
			// the corresponding nodes
			const XmlAttributeObject *typeAttr1;
			m_xmlParser->findXmlAttributeObject(*end1, type, typeAttr1);
			const XmlAttributeObject *typeAttr2;
			m_xmlParser->findXmlAttributeObject(*end2, type, typeAttr2);
			int nodeId1 = typeAttr1->m_pAttributeValue->info();
			int nodeId2 = typeAttr2->m_pAttributeValue->info();

			// Create an edge for the graph
			HashElement<int, NodeElement*> *node1HE = m_idToNode.lookup(nodeId1);
			HashElement<int, NodeElement*> *node2HE = m_idToNode.lookup(nodeId2);

			// Both nodes were found
			if (node1HE && node2HE){
				NodeElement *node1 = node1HE->info();
				NodeElement *node2 = node2HE->info();
				EdgeElement *edge = modelGraph.newEdge(node1, node2);
				modelGraph.type(edge) = Graph::association;

				// Insert edge id and edge element into hashing table
				m_idToEdge.fastInsert(edgeId, edge);
			}

			// If condition above does not hold: Error!
			// At least one node is not contained in the node hashtable
			// One reason could be that we have an association between at least
			// one element other than class or interface

			// Next association
			m_xmlParser->findBrotherXmlTagObject(*associationSon, umlAssociation, associationSon);

		} // while (associationSon != 0)

		return true;

	} // traversePackagesAndInsertAssociationEdges

	//
	// t r a v e r s e P a c k a g e s A n d I n s e r t G e n e r a l i z a t i o n E d g e s
	//
	bool UmlToGraphConverter::traversePackagesAndInsertGeneralizationEdges(
		const XmlTagObject &currentRootTag,
		UmlModelGraph &modelGraph)
	{
		// TODO: The generalization tags can also occur inside interface classifiers (in Java)
		//       Currently we only consider classes.

		// Identify contained packages (<UML:Package>)
		const XmlTagObject *packageSon;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlPackage, packageSon);
		while (packageSon != nullptr){

			// Find son umlNamespaceOwnedElement
			// if nonexistent then continue
			const XmlTagObject *newRootTag;
			m_xmlParser->findSonXmlTagObject(*packageSon, umlNamespaceOwnedElement, newRootTag);
			if (newRootTag != nullptr){

				// Call this function recursively
				if (!traversePackagesAndInsertGeneralizationEdges(*newRootTag, modelGraph))
				{
					return false;
				}

			}

			// Next package
			m_xmlParser->findBrotherXmlTagObject(*packageSon, umlPackage, packageSon);

		} // while

		// Find all classes (<UML:Class>)
		const XmlTagObject *classSon;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlClass, classSon);
		while (classSon != nullptr){

			Array<int> path(2);
			path[0] = umlNamespaceOwnedElement;
			path[1] = umlGeneralization;
			const XmlTagObject *generalizationTag = nullptr;

			// Found a <UML:Generalization> tag
			if (m_xmlParser->traversePath(*classSon, path, generalizationTag)){

				// Find out the reference number of this edge
				const XmlAttributeObject *edgeIdAttr = nullptr;
				m_xmlParser->findXmlAttributeObject(*generalizationTag, xmiId, edgeIdAttr);
				int edgeId = edgeIdAttr->m_pAttributeValue->info();

				// Find child and parent attributes
				const XmlAttributeObject *childAttr = nullptr;
				m_xmlParser->findXmlAttributeObject(*generalizationTag, child, childAttr);
				const XmlAttributeObject *parentAttr = nullptr;
				m_xmlParser->findXmlAttributeObject(*generalizationTag, parent, parentAttr);

				// Something wrong
				if (!childAttr || !parentAttr){

					GraphIO::logger.lout(Logger::LL_MINOR) << "Current dependency tag does not contain both attributes child and parent." << endl;

					// Next class
					m_xmlParser->findBrotherXmlTagObject(*classSon, umlClass, classSon);
					continue;
				}

				// Get ids and nodes
				int childId = childAttr->m_pAttributeValue->info();
				int parentId = parentAttr->m_pAttributeValue->info();

				// Get hash elements
				HashElement<int, NodeElement*> *childNodeHE = m_idToNode.lookup(childId);
				HashElement<int, NodeElement*> *parentNodeHE = m_idToNode.lookup(parentId);

				// Create an edge for the graph
				if (childNodeHE && parentNodeHE){

					NodeElement *childNode  = childNodeHE->info();
					NodeElement *parentNode = parentNodeHE->info();

					EdgeElement *edge = modelGraph.newEdge(childNode, parentNode);
					modelGraph.type(edge) = Graph::generalization;

					// Insert edge id and edge element into hashing table
					m_idToEdge.fastInsert(edgeId, edge);
				}
				// If condition above does not hold: Error!
				// At least one node is not contained in the node hashtable

			} // Found generalization tag

			// Next class
			m_xmlParser->findBrotherXmlTagObject(*classSon, umlClass, classSon);

		} // while (classSon != 0)

		return true;

	} // traversePackagesAndInsertGeneralizationEdges


	//
	// i n s e r t D e p e n d e n c y E d g e s
	//
	bool UmlToGraphConverter::insertDependencyEdges(const XmlTagObject &currentRootTag,
														UmlModelGraph &modelGraph)
	{
		// Find first dependency tag (<UML:Dependency>)
		const XmlTagObject *currentDependencyTag = nullptr;
		m_xmlParser->findSonXmlTagObject(currentRootTag, umlDependency, currentDependencyTag);

		// Find all dependencys
		while (currentDependencyTag != nullptr){

			// Find out the reference number of this edge
			const XmlAttributeObject *edgeIdAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*currentDependencyTag, xmiId, edgeIdAttr);
			int edgeId = edgeIdAttr->m_pAttributeValue->info();

			// Find client and supplier attributes
			const XmlAttributeObject *clientAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*currentDependencyTag, client, clientAttr);
			const XmlAttributeObject *supplierAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*currentDependencyTag, supplier, supplierAttr);

			// Something wrong
			if (!clientAttr || !supplierAttr){

				GraphIO::logger.lout(Logger::LL_MINOR) << "Current dependency tag does not contain both attributes client and supplier." << endl;

				// Next dependency
				m_xmlParser->findBrotherXmlTagObject(*currentDependencyTag, umlDependency, currentDependencyTag);
				continue;
			}

			// Get ids
			int clientId = clientAttr->m_pAttributeValue->info();
			int supplierId = supplierAttr->m_pAttributeValue->info();

			// Get Hashelements
			HashElement<int, NodeElement*> *clientNodeHE = m_idToNode.lookup(clientId);
			HashElement<int, NodeElement*> *supplierNodeHE = m_idToNode.lookup(supplierId);

			// Create an edge for the graph
			if (clientNodeHE && supplierNodeHE){

				NodeElement *clientNode   = clientNodeHE->info();
				NodeElement *supplierNode = supplierNodeHE->info();

				EdgeElement *edge = modelGraph.newEdge(clientNode, supplierNode);
				modelGraph.type(edge) = Graph::dependency;

				// Insert edge id and edge element into hashing table
				m_idToEdge.fastInsert(edgeId, edge);
			}
			// If condition above does not hold: Error!
			// At least one node is not contained in the node hashtable

			// Next dependency
			m_xmlParser->findBrotherXmlTagObject(*currentDependencyTag, umlDependency, currentDependencyTag);

		} // while (currentDependecyTag != 0)

		return true;

	} // insertDependencyEdges


	//
	// s t r i n g T o D o u b l e A r r a y
	//
	// Extracts the single values of string str with format
	// "x, y, width, height," and puts them into doubleArray
	static void stringToDoubleArray(const string &str, Array<double> &doubleArray)
	{
		size_t strIndex = 0;
		char tempString[20];

		for (int i = 0; i < 4; i++){

			int tempStringIndex = 0;

			// Skip whitespace
			while (isspace(str[strIndex])){
				++strIndex;
			}

			// Copy characters of double value
			// values are separated by comma
			while (str[strIndex] != ','){

				tempString[tempStringIndex] = str[strIndex];
				++tempStringIndex;
				++strIndex;
			}

			// Skip over ','
			++strIndex;

			// Terminate string
			tempString[tempStringIndex] = '\0';

			// Put double value into array
			doubleArray[i] = atof(tempString);

		} // for

	} // stringToDoubleArray


	//
	// c r e a t e D i a g r a m G r a p h s
	//
	bool UmlToGraphConverter::createDiagramGraphs()
	{
		// We want to create a diagram graph for each subtree <UML:Diagram> found
		// in the parse tree.
		//
		// Currently we are only interested in class diagrams.

		// Model graph must exist!
		OGDF_ASSERT(m_modelGraph != nullptr);

		// Message
		//cout << "Creating diagram graph(s)..." << endl;

		// Check root element (must be <XMI>)
		if (m_xmlParser->getRootTag().m_pTagName->info() != xmi){
			GraphIO::logger.lout() << "Root tag is not <XMI>" << endl;
			return false;
		}

		// Find the first <UML:Diagram> tag starting at <XMI>
		Array<int> path(2);
		path[0] = xmiContent;
		path[1] = umlDiagram;
		const XmlTagObject *currentDiagramTag = nullptr;
		m_xmlParser->traversePath(m_xmlParser->getRootTag(), path, currentDiagramTag);

		// Traverse diagrams
		while (currentDiagramTag != nullptr){

			// Find out name of the diagram
			const XmlAttributeObject *nameAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*currentDiagramTag, name, nameAttr);
			string diagramName("");
			if (nameAttr != nullptr){
				diagramName = nameAttr->m_pAttributeValue->key();
			}

			// Find out type of the diagram
			const XmlAttributeObject *diagramTypeAttr = nullptr;
			m_xmlParser->findXmlAttributeObject(*currentDiagramTag, diagramType, diagramTypeAttr);

			// No diagramTypeAttribute found --> we continue with the next diagram
			if (diagramTypeAttr == nullptr){

				// Next diagram
				m_xmlParser->findBrotherXmlTagObject(*currentDiagramTag, umlDiagram, currentDiagramTag);
				continue;
			}

			// Check which type of diagram we have
			UmlDiagramGraph::UmlDiagramType diagramType;
			switch (diagramTypeAttr->m_pAttributeValue->info()){

			case (classDiagram) :
				diagramType = UmlDiagramGraph::classDiagram;
				break;
			case (moduleDiagram) :
				diagramType = UmlDiagramGraph::moduleDiagram;
				break;
			default:
				diagramType = UmlDiagramGraph::unknownDiagram;
				break;

			} // switch

			// Currently we only allow class diagrams; in all other cases
			// we continue with the next diagram
			if (diagramType != UmlDiagramGraph::classDiagram){

				// Next diagram
				m_xmlParser->findBrotherXmlTagObject(*currentDiagramTag, umlDiagram, currentDiagramTag);
				continue;
			}

			// Create a new diagram graph and add it to the list of diagram graphs
			UmlDiagramGraph *diagramGraph =
				new UmlDiagramGraph(*m_modelGraph,
										diagramType,
										diagramName);
			m_diagramGraphs.pushBack(diagramGraph);


			// First pass the <UML:Diagram.element> tag
			const XmlTagObject *rootDiagramElementTag = nullptr;
			m_xmlParser->findSonXmlTagObject(*currentDiagramTag, rootUmlDiagramElement, rootDiagramElementTag);

			// No such tag found --> we continue with the next diagram
			if (rootDiagramElementTag == nullptr){

				// Next diagram
				m_xmlParser->findBrotherXmlTagObject(*currentDiagramTag, umlDiagram, currentDiagramTag);
				continue;
			}

			// Now investigate the diagram elements
			const XmlTagObject *currentDiagramElementTag = nullptr;
			m_xmlParser->findSonXmlTagObject(*rootDiagramElementTag, umlDiagramElement, currentDiagramElementTag);

			// Traverse all diagram elements (<UML:DiagramElement>)
			while (currentDiagramElementTag != nullptr){

				// We have to investigate the subject attribute which contains the
				// reference number of the represented element; then we can check if
				// a node or edge with this reference exists
				const XmlAttributeObject *subjectAttr = nullptr;
				m_xmlParser->findXmlAttributeObject(*currentDiagramElementTag, subject, subjectAttr);

				// Not found --> continue with the next diagram element
				if (subjectAttr == nullptr){

					// Next diagram element
					m_xmlParser->findBrotherXmlTagObject(*currentDiagramElementTag, umlDiagramElement, currentDiagramElementTag);

					continue;
				}

				// Check wether node or edge with this reference does exist
				int elementId = subjectAttr->m_pAttributeValue->info();

				// Node exists for that reference
				if (m_idToNode.lookup(elementId) != nullptr){

					// Get hash element
					HashElement<int,NodeElement*> *nodeHashElement = m_idToNode.lookup(elementId);

					// Get node element
					NodeElement* geometricNode = nodeHashElement->info();

					// Extract geometric information
					const XmlAttributeObject *geometryAttr = nullptr;
					m_xmlParser->findXmlAttributeObject(*currentDiagramElementTag, geometry, geometryAttr);

					// Not found
					if (geometryAttr == nullptr){
						break;
					}

					// Get double values of geometry
					Array<double> geometryArray(4);
					stringToDoubleArray(geometryAttr->m_pAttributeValue->key(), geometryArray);

					// Add node to diagram graph
					diagramGraph->addNodeWithGeometry(
						geometricNode,
						geometryArray[0],
						geometryArray[1],
						geometryArray[2],
						geometryArray[3]);

				} // Node exists
				// Node does not exist
				else{

					// Edge exists for that reference
					if (m_idToEdge.lookup(elementId) != nullptr){

						// Get hash element
						HashElement<int,EdgeElement*> *edgeHashElement = m_idToEdge.lookup(elementId);

						// Get node element
						EdgeElement* geometricEdge = edgeHashElement->info();

						// Add edge to diagram graph
						diagramGraph->addEdge(geometricEdge);

					} // Edge exists

				} // else

				// Next diagram element
				m_xmlParser->findBrotherXmlTagObject(*currentDiagramElementTag,
													 umlDiagramElement,
													 currentDiagramElementTag);

			} // while (currentDiagramElementTag != 0)

			// Next diagram
			m_xmlParser->findBrotherXmlTagObject(*currentDiagramTag, umlDiagram, currentDiagramTag);


		} // while (currentDiagramTag != 0)

		return true;

	} // createDiagramGraphs

	//
	// c r e a t e D i a g r a m G r a p h s I n U M L G r a p h F o r m a t
	//
	bool UmlToGraphConverter::createDiagramGraphsInUMLGraphFormat(SList<UMLGraph*> &diagramGraphsInUMLGraphFormat)
	{
		// We want to create an instance of UMLGraph for each instance of UmlDiagramGraph
		// contained in the given list. Implicitly we have to create also an instance of class Graph
		// for each UMLGraph.
		// We maintain a hash list for mapping the nodes and edges of the model graph to the
		// new nodes and edges of the graph created for the diagram.
		// We use as key the unique index of the node resp. edge.

		// Message
		//cout << "Creating diagram graph(s) in UMLGraph format..." << endl;

		// Traverse list of diagram graphs
		for (UmlDiagramGraph *diagramGraph : m_diagramGraphs)
		{
			// Mapping from the index of the existing node to the new nodeElement
			Hashing<int, NodeElement*> indexToNewNode;

			// Mapping from the index of the existing edge to the new edgeElement
			Hashing<int, EdgeElement*> indexToNewEdge;

			// Create instance of class graph
			Graph *graph = new Graph();

			// Traverse list of nodes contained in the diagram
			const SList<NodeElement*> &diagramNodes = diagramGraph->getNodes();
			for (node n : diagramNodes)
			{
				// Create a new "pendant" node for the existing	node
				NodeElement *newNode = graph->newNode();

				// Insert mapping from index of the existing node to the pendant node
				// into hashtable
				indexToNewNode.fastInsert(n->index(), newNode);

			} // Traverse list of nodes contained in the diagram

			// Traverse list of edges contained in the diagram
			const SList<EdgeElement*> &diagramEdges = diagramGraph->getEdges();
			for (edge e : diagramEdges)
			{
				// Find out source and target of the edge
				NodeElement *source = e->source();
				NodeElement *target = e->target();

				// Find pendant nodes
				HashElement<int, NodeElement*> *sourceHashElement =
					indexToNewNode.lookup(source->index());
				HashElement<int, NodeElement*> *targetHashElement =
					indexToNewNode.lookup(target->index());
				NodeElement *pendantSource = sourceHashElement->info();
				NodeElement *pendantTarget = targetHashElement->info();

				// Insert new edge between pendant nodes
				EdgeElement *newEdge = graph->newEdge(pendantSource, pendantTarget);

				// Insert mapping from index of the existing edgeto the pendant edge
				// into hashtable
				indexToNewEdge.fastInsert(e->index(), newEdge);

			} // Traverse list of edges contained in the diagram

			// Create instance of class UMLGraph
			UMLGraph *umlGraph = new UMLGraph(*graph, GraphAttributes::nodeLabel);

			// Now we want to add the geometry information and the node label
			const SList<double> xList = diagramGraph->getX();
			const SList<double> yList = diagramGraph->getY();
			const SList<double> wList = diagramGraph->getWidth();
			const SList<double> hList = diagramGraph->getHeight();
			SListConstIterator<double> xIt, yIt, wIt, hIt;

			// Traverse node list and geometry lists synchronously
			xIt = xList.begin();
			yIt = yList.begin();
			wIt = wList.begin();
			hIt = hList.begin();
			for(node n : diagramNodes)
			{
				// Get pendant node
				HashElement<int, NodeElement*> *nodeHashElement =
					indexToNewNode.lookup(n->index());
				NodeElement *pendantNode = nodeHashElement->info();

				// Insert geometry information into umlGraph
				umlGraph->x(pendantNode) = *xIt;
				umlGraph->y(pendantNode) = *yIt;
				umlGraph->width(pendantNode) = *wIt;
				umlGraph->height(pendantNode) = *hIt;

				// Insert label
				string &label = umlGraph->label(pendantNode);
				label = m_modelGraph->getNodeLabel(n);

				// Next iteration
				++xIt;
				++yIt;
				++wIt;
				++hIt;

			} // Traverse node list and geometry lists synchronously

			// Traverse list of edges contained in the diagram
			for (edge e : diagramEdges)
			{
				// Find pendant edge
				HashElement<int, EdgeElement*> *edgeHashElement =
					indexToNewEdge.lookup(e->index());
				EdgeElement *pendantEdge = edgeHashElement->info();

				// Insert type information into umlGraph
				umlGraph->type(pendantEdge) = m_modelGraph->type(e);

			} // Traverse edge list and insert type information for the new edges.

			// Add new umlGraph to list
			diagramGraphsInUMLGraphFormat.pushBack(umlGraph);

		} // Traverse list of diagram graphs

		return true;

	} // createDiagramGraphsInUMLGraphFormat

} // namespace ogdf
