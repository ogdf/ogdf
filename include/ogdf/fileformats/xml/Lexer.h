/** \file
 * \brief Declarations for simple XML lexer.
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


#ifndef OGDF_XML_LEXER_H
#define OGDF_XML_LEXER_H


#include <iostream>
#include <array>
#include <string>


namespace ogdf {


namespace xml {


enum class token {
	chevron_left,
	chevron_right,
	questionMark,
	assignment,
	slash,
	string,
	identifier,
	eof
};


class Lexer {
public:
	Lexer();
	Lexer(std::istream &is);

	void setInput(std::istream &is);

	const std::string &currentValue() const {
		return m_value;
	}

	token currentToken() const {
		return m_token;
	}

	token nextToken();
	void consumeText();

	std::size_t row() const {
		return m_row;
	}

	std::size_t column() const {
		return m_column;
	}

private:
	std::istream *m_input;
	std::size_t m_row;
	std::size_t m_column;

	token m_token;
	std::string m_value;

	std::array<char, 4> m_buffer;
	std::size_t m_buffer_begin;
	std::size_t m_buffer_size;

	std::size_t buffer_size() const;
	bool buffer_empty() const;
	char buffer_peek() const;
	char buffer_at(std::size_t i) const;
	void buffer_put(int c);
	char buffer_pop();

	void skipTrash();
	bool skipWhitespace();
	bool skipComments();

	void consumeString();
	void consumeIdentifier();
};


inline std::size_t Lexer::buffer_size() const
{
	return m_buffer_size;
}


inline bool Lexer::buffer_empty() const
{
	return buffer_size() == 0;
}


inline char Lexer::buffer_peek() const
{
	return m_buffer[m_buffer_begin];
}


inline char Lexer::buffer_at(std::size_t ind) const
{
	return m_buffer[(m_buffer_begin + ind) % m_buffer.size()];
}


inline void Lexer::buffer_put(int c)
{
	m_buffer[(m_buffer_begin + buffer_size()) % m_buffer.size()] = static_cast<char>(c);
	m_buffer_size++;
}


inline char Lexer::buffer_pop()
{
	char result = buffer_peek();
	m_buffer_begin = (m_buffer_begin + 1) % (m_buffer.size());
	m_buffer_size--;

	// Adjust cursor position.
	m_column++;
	if (result == '\n') {
		m_row++;
		m_column = 0;
	}

	return result;
}


}


}


#endif
