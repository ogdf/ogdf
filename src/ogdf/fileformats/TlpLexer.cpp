/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of TLP file format lexer class.
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

#include <ogdf/fileformats/TlpLexer.h>


namespace ogdf {

namespace tlp {


std::ostream &operator <<(std::ostream &os, const Token &token)
{
	switch(token.type) {
	case Token::tok_leftParen:
		os << "tok_(";
		break;
	case Token::tok_rightParen:
		os << "tok_)";
		break;
	case Token::tok_identifier:
		os << "tok_id(" << *(token.value) << ")";
		break;
	case Token::tok_string:
		os << "tok_str(\"" << *(token.value) << "\")";
		break;
	}

	return os;
}


Token::Token(const Type &type, size_t line, size_t column)
: type(type), line(line), column(column)
{
	if(type == tok_identifier || type == tok_string) {
		value = new std::string;
	} else {
		value = NULL;
	}
}


bool Lexer::isIdentifier(char c)
{
	return isalnum(c) || c == '_' || c == '.' || c == '-';
}


bool Lexer::fetchBuffer()
{
	if(!std::getline(m_istream, m_buffer)) {
		return false;
	}

	m_begin = m_buffer.begin();
	m_end = m_buffer.end();

	m_line++;
	return true;
}


void Lexer::cleanValues()
{
	for(std::vector<Token>::iterator it = m_tokens.begin();
	    it != m_tokens.end();
	    it++)
	{
		delete it->value;
	}
}


Lexer::Lexer(std::istream &is) : m_istream(is)
{
}


Lexer::~Lexer()
{
	cleanValues();
}


bool Lexer::tokenize()
{
	cleanValues();
	m_tokens.clear();

	m_line = 0;
	while(fetchBuffer()) {
		if(!tokenizeLine()) {
			return false;
		}
	}

	return true;
}


bool Lexer::tokenizeLine()
{
	while(m_begin != m_end && isspace(*m_begin)) {
		m_begin++;
	}

	// We got an end of a line or a comment.
	if(m_begin == m_end || *m_begin == ';') {
		return true;
	}

	if(*m_begin == '(') {
		m_tokens.push_back(Token(Token::tok_leftParen, line(), column()));
		m_begin++;
		return tokenizeLine();
	}

	if(*m_begin == ')') {
		m_tokens.push_back(Token(Token::tok_rightParen, line(), column()));
		m_begin++;
		return tokenizeLine();
	}

	if(*m_begin == '"') {
		return tokenizeString() && tokenizeLine();
	}

	if(isIdentifier(*m_begin)) {
		return tokenizeIdentifier() && tokenizeLine();
	}

	std::cerr << "ERROR: Unexpected character \"" << *m_begin << "\" at ("
	          << line() << ", " << column() << ").\n";
	return false;
}


bool Lexer::tokenizeString()
{
	m_begin++;

	Token token(Token::tok_string, line(), column());

	for(;;) {
		// Check whether we need to refill the buffer.
		if(m_begin == m_end && !fetchBuffer()) {
			std::cerr << "ERROR: End of input while parsing a string at ("
			          << token.line << ", "
			          << token.column << ").\n";
			return false;
		}

		if(m_begin == m_end) {
			continue;
		}

		// Check whether we got a end of a string.
		if(*m_begin == '"') {
			m_tokens.push_back(token);
			m_begin++;
			return true;
		}

		// If the to above failed, we just put new character.
		*(token.value) += *m_begin;
		m_begin++;
	}

	return true;
}


bool Lexer::tokenizeIdentifier()
{
	Token token(Token::tok_identifier, line(), column());

	while(m_begin != m_end && isIdentifier(*m_begin)) {
		*(token.value) += *m_begin;
		m_begin++;
	}

	m_tokens.push_back(token);
	return true;
}


} // end namespace tlp

} // end namespace ogdf

