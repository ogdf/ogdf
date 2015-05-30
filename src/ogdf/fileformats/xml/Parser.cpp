/** \file
 * \brief Implementation of simple XML parser.
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

#include <ogdf/fileformats/xml/Parser.h>

#include <sstream>


namespace ogdf {


namespace xml {


ParseException::ParseException(
	const std::string &what,
	std::size_t row, std::size_t column)
{
	std::ostringstream message;
	message << row << ":" << column << ": " << what;
	m_what = message.str();
}


ParseException::ParseException(
	const char *what,
	std::size_t row, std::size_t column)
: ParseException(std::string(what), row, column)
{
}


unique_ptr<Element> Parser::parse(std::istream &is)
{
	m_lexer.setInput(is);

	Piece piece = std::move(nextPiece());
	if (piece.m_type == Piece::type::tag_header) {
		// TODO: Extract some useful information from XML header.
		piece = std::move(nextPiece());
	}

	switch (piece.m_type) {
	case Piece::type::tag_open:
	case Piece::type::tag_empty:
		return std::move(parseElement(piece));
	default:
		fail("expected root tag");
		return std::move(unique_ptr<Element>(nullptr));
	}
}


unique_ptr<Element> Parser::parseElement(Piece &piece)
{
	unique_ptr<Element> elem(new Element());
	elem->m_name = std::move(piece.m_value);
	elem->m_attrs = std::move(*piece.m_attrs);

	if (piece.m_type == Piece::type::tag_empty) {
		return std::move(elem);
	}

	for (;;) {
		Piece next = nextPiece();
		switch (next.m_type) {
		case Piece::type::text:
			elem->m_text.append(std::move(next.m_value));
			break;
		case Piece::type::tag_open:
		case Piece::type::tag_empty:
			elem->m_children.push_back(std::move(parseElement(next)));
			break;
		case Piece::type::tag_close:
			if (next.m_value != elem->m_name) {
				fail("expected tag closing \"" + elem->m_name + "\" " +
					 "but found \"" + next.m_value);
			}
			return std::move(elem);
		case Piece::type::tag_header:
			/*
			 * Header should happen only as the first statement so I guess it
			 * should fail here. But simply ignoring it means no charm, right?
			 */
			break;
		}
	}
}


Parser::Piece Parser::nextPiece()
{
	Piece piece;

	m_lexer.consumeText();
	if (!isBlank(m_lexer.currentValue())) {
		piece.m_type = Piece::type::text;
		piece.m_value = m_lexer.currentValue();
		return std::move(piece);
	}

	m_lexer.nextToken();
	expect(token::chevron_left, "expected opening tag bracket ('<')");

	m_lexer.nextToken();
	if (m_lexer.currentToken() == token::slash) {
		m_lexer.nextToken();
		expect(token::identifier, "expected close tag name");

		piece.m_type = Piece::type::tag_close;
		piece.m_value = m_lexer.currentValue();

		m_lexer.nextToken();
		expect(token::chevron_right, "exected closing tag bracket ('>')");
		return std::move(piece);
	}
	if (m_lexer.currentToken() == token::questionMark) {
		m_lexer.nextToken();

		piece.m_type = Piece::type::tag_header;
	}

	expect(token::identifier, "expected tag name");

	piece.m_type = Piece::type::tag_open;
	piece.m_value = m_lexer.currentValue();
	piece.m_attrs.reset(new Element::Attributes());

	for (;;) {
		m_lexer.nextToken();
		if (m_lexer.currentToken() == token::slash) {
			m_lexer.nextToken();
			expect(token::chevron_right, "expected closing tag bracket ('>')");

			piece.m_type = Piece::type::tag_empty;
		}
		if (m_lexer.currentToken() == token::questionMark) {
			m_lexer.nextToken();
			expect(token::chevron_right, "expected closing tag bracket ('>')");

			piece.m_type = Piece::type::tag_header;
		}
		if (m_lexer.currentToken() == token::chevron_right) {
			return std::move(piece);
		}

		expect(token::identifier, "expected attribute identifier");
		std::string key = m_lexer.currentValue();

		m_lexer.nextToken();
		expect(token::assignment, "expected assignment operator ('=')");

		m_lexer.nextToken();
		expect(token::string, "exected string value of attribute");
		std::string value = m_lexer.currentValue();

		(*piece.m_attrs)[key] = value;
	}
}


}


}
