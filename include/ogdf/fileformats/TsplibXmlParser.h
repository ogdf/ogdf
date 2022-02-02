/** \file
 * \brief Parser for TSPLIB instances in XML format
 *
 * \author Finn Stutzenstein
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

#include <ogdf/fileformats/GraphIO.h>

#include <ogdf/basic/HashArray.h>
#include <ogdf/lib/pugixml/pugixml.h>

#include <sstream>
#include <unordered_map>

namespace ogdf {

/**
 * Parses tsplib files in xml format.
 *
 * See http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/XML-TSPLIB/Description.pdf
 *
 * Notes:
 * - Node indices will be set to the indices of the instance, so all nodes
 *   will have indices from 0 to n-1.
 * - Selfloops will be skipped
 */
class TsplibXmlParser {
private:
	pugi::xml_document m_xml; // hold as a class member, so it gets desctructed at the same time as the parser.
	pugi::xml_node m_graphTag;

	bool m_hasError;

	// Called on constructing, to initially load the file.
	// Basic errors are catched, like missing essential tags.
	// Parsing is done on read, which might also raise errors.
	bool load(std::istream& in, std::string& error);

	// Unified read with all the logic.
	bool read(Graph& G, GraphAttributes* GA);
public:
	explicit TsplibXmlParser(std::istream& in);
	~TsplibXmlParser() = default;

	bool read(Graph& G);
	bool read(Graph& G, GraphAttributes& GA);
};

}
