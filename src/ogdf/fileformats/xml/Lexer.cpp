/** \file
 * \brief Implementation of simple XML lexer.
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

#include <ogdf/fileformats/xml/Lexer.h>


namespace ogdf {


namespace xml {


Lexer::Lexer()
: m_input(nullptr)
{
}


Lexer::Lexer(std::istream &is)
{
	setInput(is);
}


void Lexer::setInput(std::istream &is)
{
	m_input = &is;

	// Input changes, so we need to reset our position...
	m_row = 1;
	m_column = 0;

	// ... and reset the buffer.
	m_buffer_begin = 0;
	m_buffer_size = 0;
}


xml::token Lexer::nextToken()
{
	while (skipWhitespace() || skipComments()) {
	}

	if (buffer_empty()) {
		buffer_put(m_input->get());
	}

	switch (buffer_peek()) {
	case '<':
		m_token = token::chevron_left;
		buffer_pop();
		break;
	case '>':
		m_token = token::chevron_right;
		buffer_pop();
		break;
	case '?':
		m_token = token::questionMark;
		buffer_pop();
		break;
	case '=':
		m_token = token::assignment;
		buffer_pop();
		break;
	case '/':
		m_token = token::slash;
		buffer_pop();
		break;
	case '"': case '\'':
		consumeString();
		break;
	case EOF:
		/*
		 * TODO: Replace EOF with std::char_traits<char>::eof() once VS
		 * will learn what "constexpr" is.
		 */
		m_token = token::eof;
		buffer_pop();
		break;
	default:
		consumeIdentifier();
		break;
	}

	return m_token;
}


bool Lexer::skipWhitespace()
{
	bool consumed = false;

	for (;;) {
		if (buffer_empty()) {
			buffer_put(m_input->get());
		}
		char current = buffer_peek();

		// TODO: Replace it with std::isspace in the future.
		if (!isspace(current)) {
			break;
		}

		buffer_pop();
		consumed = true;
	}

	return consumed;
}


bool Lexer::skipComments()
{
	// We need to ensure that we have enough characters in the buffer.
	switch (buffer_size()) {
	case 0:
		buffer_put(m_input->get());
	case 1:
		buffer_put(m_input->get());
	case 2:
		buffer_put(m_input->get());
	case 3:
		buffer_put(m_input->get());
	}

	// We have a comment iff it begins with '<!--' sequence.
	if (!(buffer_at(0) == '<' &&
	      buffer_at(1) == '!' &&
	      buffer_at(2) == '-' &&
	      buffer_at(3) == '-'))
	{
		return false;
	}

	buffer_pop();
	buffer_pop();
	buffer_pop();
	buffer_pop();

	for (;;) {
		// TODO: Handle unclosed comments.

		// As above, we enusre that we have enough characters available.
		switch (buffer_size()) {
		case 0:
			buffer_put(m_input->get());
		case 1:
			buffer_put(m_input->get());
		case 2:
			buffer_put(m_input->get());
		}

		// The comment ends only with the '-->' sequence.
		if (buffer_at(0) == '-' &&
		    buffer_at(1) == '-' &&
		    buffer_at(2) == '>')
		{
			buffer_pop();
			buffer_pop();
			buffer_pop();
			break;
		}

		buffer_pop();
	}

	return true;
}


void Lexer::consumeText()
{
	m_value = "";

	for (;;) {
		while (skipComments()) {
		}

		if (buffer_empty()) {
			buffer_put(m_input->get());
		}
		char current = buffer_peek();

		// TODO: Escape '&lt;' and '&gt;' as '<' and '>';
		if (current == '<' || current == '>' ||
		    current == std::char_traits<char>::eof())
		{
			break;
		}

		m_value.push_back(buffer_pop());
	}
}


void Lexer::consumeString()
{
	m_value = "";

	char delim = buffer_pop();
	for (;;) {
		if (buffer_empty()) {
			buffer_put(m_input->get());
		}
		char current = buffer_peek();

		if (current == delim) {
			buffer_pop();
			break;
		}

		if (current == std::char_traits<char>::eof()) {
			break;
		}

		m_value.push_back(buffer_pop());
	}

	m_token = token::string;
}


void Lexer::consumeIdentifier()
{
	m_value = "";

	for (;;) {
		if (buffer_empty()) {
			buffer_put(m_input->get());
		}
		char current = buffer_peek();

		/*
		 * The XML spec allows much narrower set of valid identifiers, but
		 * I see not charm in accepting some of the malformed identifiers
		 * since the parser is not very strict about standard too.
		 */
		// TODO: Replace it with std::isspace in the future.
		if (isspace(current) ||
		    current == '<' || current == '>' ||
		    current == '=' || current == '/' || current == '?' ||
		    current == '"' || current == '\'' ||
		    current == std::char_traits<char>::eof())
		{
			break;
		}

		m_value.push_back(buffer_pop());
	}

	m_token = token::identifier;
}


}


}
