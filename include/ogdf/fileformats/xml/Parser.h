/** \file
 * \brief Declarations for simple XML parser.
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


#ifndef OGDF_XML_PARSER_H
#define OGDF_XML_PARSER_H


#include <algorithm>
#include <exception>
#include <memory>
#include <vector>
#include <unordered_map>
using std::unique_ptr;

#include <ogdf/fileformats/xml/Lexer.h>


namespace ogdf {


namespace xml {


class ParseException : public std::exception {
public:
	explicit ParseException(
		const std::string &what,
		std::size_t row, std::size_t column);
	explicit ParseException(
		const char *what,
		std::size_t row, std::size_t column);

	virtual const char *what() const throw() override {
		return m_what.c_str();
	}

private:
	std::string m_what;
};


class Element {
public:
	using Children = std::vector<unique_ptr<Element>>;
	using Attributes = std::unordered_map<std::string, std::string>;

	const std::string &name() const {
		return m_name;
	}

	const std::string &text() const {
		return m_text;
	}

	const Children &children() const {
		return m_children;
	}

	const Attributes &attributes() const {
		return m_attrs;
	}

private:
	std::string m_name;
	std::string m_text;
	Children m_children;
	Attributes m_attrs;

	friend class Parser;
};


class Parser {
public:
	unique_ptr<Element> parse(std::istream &is);

private:
	Lexer m_lexer;

	class Piece {
	public:
		enum class type {
			text,
			tag_open,
			tag_close,
			tag_empty,
			tag_header
		};

		type m_type;
		std::string m_value;
		std::unique_ptr<Element::Attributes> m_attrs;
	};

	Piece nextPiece();
	std::unique_ptr<Element> parseElement(Piece &piece);

	void fail(const char *message);
	void fail(const std::string &message);
	void expect(token expected, const char *message);

	static bool isBlank(const std::string &str);
};


inline void Parser::fail(const char *message)
{
	throw ParseException(message, m_lexer.row(), m_lexer.column());
}


inline void Parser::fail(const std::string &message)
{
	fail(message.c_str());
}


inline void Parser::expect(token expected, const char *message)
{
	if (m_lexer.currentToken() == expected) {
		return;
	}

	fail(message);
}


inline bool Parser::isBlank(const std::string &str)
{
	return std::all_of(str.begin(), str.end(), [](char c) {
		// TODO: Replace it with std::isspace once VS will stop being retarded.
		return isspace(c);
	});
}


}


}


#endif
