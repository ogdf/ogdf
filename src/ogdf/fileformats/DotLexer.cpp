/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements DOT format Lexer class.
 *
 * \author Łukasz Hanuszczak
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).

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

#include <ogdf/fileformats/DotLexer.h>

namespace ogdf {

namespace dot {


Token::Token(
	size_t row, size_t column,
	std::string *value)
: row(row), column(column), value(value)
{
}


std::string Token::toString(const Type &type)
{
	switch(type) {
	case assignment: return "=";
	case colon: return ":";
	case semicolon: return ";";
	case comma: return ",";
	case edgeOpDirected: return "->";
	case edgeOpUndirected: return "--";
	case leftBracket: return "[";
	case rightBracket: return "]";
	case leftBrace: return "{";
	case rightBrace: return "}";
	case graph: return "graph";
	case digraph: return "digraph";
	case subgraph: return "subgraph";
	case node: return "node";
	case edge: return "edge";
	case strict: return "strict";
	case identifier: return "identifier";
	default: return "unknown";
	}
}


Lexer::Lexer(std::istream &input) : m_input(input)
{
}


Lexer::~Lexer()
{
	for(std::vector<Token>::iterator it = m_tokens.begin();
	    it != m_tokens.end();
	    it++)
	{
		delete it->value;
	}
}


const std::vector<Token> &Lexer::tokens() const
{
	return m_tokens;
}


bool Lexer::tokenize()
{
	m_row = 0;
	while(m_input.good()) {
		if(!tokenizeLine()) {
			return false;
		}
	}

	return true;
}


bool Lexer::tokenizeLine()
{
	std::getline(m_input, m_buffer);
	m_row++;

	// Handle line output from a C preprocessor (#blabla).
	if(m_buffer[0] == '#') {
		return true;
	}

	for(m_col = 0; m_col < m_buffer.size(); m_col++) {
		// Ignore whitespaces.
		if(isspace(m_buffer[m_col])) {
			continue;
		}

		// Handle single-line comments.
		if(match("//")) {
			break;
		}

		// Handle multi-line comments.
		if(match("/*")) {
			const size_t column = m_col;
			const size_t row = m_row;

			do {
				m_col++;

				// Get a new line if a current one has ended.
				if(m_col >= m_buffer.size()) {
					if(!m_input.good()) {
						std::cerr << "ERROR: Unclosed comment at"
						          << column << ", " << row;
						return false;
					}
					std::getline(m_input, m_buffer);
					m_row++;
					m_col = 0;
				}
			} while(!(m_buffer[m_col - 1] == '*' && m_buffer[m_col] == '/'));

			m_col += 2;
			continue;
		}

		Token token(m_row + 1, m_col + 1);

		if(match(Token::assignment)) {
			token.type = Token::assignment;
		} else if(match(Token::colon)) {
			token.type = Token::colon;
		} else if(match(Token::semicolon)) {
			token.type = Token::semicolon;
		} else if(match(Token::comma)) {
			token.type = Token::comma;
		} else if(match(Token::edgeOpDirected)) {
			token.type = Token::edgeOpDirected;
		} else if(match(Token::edgeOpUndirected)) {
			token.type = Token::edgeOpUndirected;
		} else if(match(Token::leftBracket)) {
			token.type = Token::leftBracket;
		} else if(match(Token::rightBracket)) {
			token.type = Token::rightBracket;
		} else if(match(Token::leftBrace)) {
			token.type = Token::leftBrace;
		} else if(match(Token::rightBrace)) {
			token.type = Token::rightBrace;
		} else if(match(Token::graph)) {
			token.type = Token::graph;
		} else if(match(Token::digraph)) {
			token.type = Token::digraph;
		} else if(match(Token::subgraph)) {
			token.type = Token::subgraph;
		} else if(match(Token::node)) {
			token.type = Token::node;
		} else if(match(Token::edge)) {
			token.type = Token::edge;
		} else if(match(Token::strict)) {
			token.type = Token::strict;
		} else if(identifier(token)) {
			token.type = Token::identifier;
		} else {
			std::cerr << "EROR: Unknown token at: "
			          << m_row << "; " << m_col
			          << "\n";
			return false;
		}

		m_tokens.push_back(token);
	}

	return true;
}


bool Lexer::match(const Token::Type &type)
{
	return match(Token::toString(type));
}


bool Lexer::match(const std::string &str)
{
	// Check whether buffer is too short to match.
	if(m_buffer.length() - m_col < str.length()) {
		return false;
	}

	for(size_t i = 0; i < str.length(); i++) {
		if(m_buffer[m_col + i] != str[i]) {
			return false;
		}
	}

	// After successful match we move the "head".
	m_col += str.length() - 1;

	return true;
}


bool Lexer::identifier(Token &token)
{
	// Check whether identifier is double-quoted string.
	if(m_buffer[m_col] == '"') {
		m_col++;
		std::stringstream ss;

		while(m_buffer[m_col] != '"' || m_buffer[m_col - 1] == '\\') {
			ss << m_buffer[m_col++];

			// Get a new line if a current one has ended.
			if(m_col >= m_buffer.size()) {
				if(!m_input.good()) {
					std::cerr << "ERROR: Unclosed string at "
					          << token.row << ", " << token.column
					          << ".\n";
					return false;
				}
				std::getline(m_input, m_buffer);
				m_row++;
				m_col = 0;
			}
		}

		token.value = new std::string(ss.str());
		return true;
	}

	// Check whether identifier is a normal C-like identifier.
	if(isalpha(m_buffer[m_col]) || m_buffer[m_col] == '_') {
		std::ostringstream ss;

		while(isalnum(m_buffer[m_col]) || m_buffer[m_col] == '_') {
			ss << m_buffer[m_col++];
		}

		m_col--;
		token.value = new std::string(ss.str());
		return true;
	}

	// Check whether identifier is a numeric literal. Quite ugly and slow but works.
	std::istringstream ss(m_buffer.c_str() + m_col);
	double temp;
	if(ss >> temp) {
		size_t length = static_cast<size_t>(ss.tellg());
		token.value = new std::string(m_buffer.substr(m_col, length));
		m_col += length;
		return true;
	}

	// TODO: HTML string identifiers.

	return false;
}


} // end namespace dot

} // end namespace ogdf
