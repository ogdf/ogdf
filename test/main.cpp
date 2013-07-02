/*
 * $Revision$
 *
 * last checkin:
 *   $Author$
 *   $Date$
 ***************************************************************/

/** \file
 * \brief Implementation of a command line based tool to run
 * tests.
 *
 * \author Christoph Schulz, Stephan Beyer
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


#include "gtest/gtest.h"
extern int regressionMain(int argc, char *const argv[]);

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2 || strlen(argv[1]) != 2 || argv[1][0] !=  '-') {
		cout
		  << "Usage: " << argv[0] << " -G [GTest options]\n"
		  << "   or: " << argv[0] << " -R [Regression test index]\n";
		return -1;
	}
	--argc;
	++argv;
	switch (argv[0][1]) {
	case 'G':
		::testing::InitGoogleTest(&argc, argv);
		return RUN_ALL_TESTS();
	case 'R':
		return regressionMain(argc, argv);
	default:
		cerr << "Error: first argument has to be either -G or -R!\n";
		return -1;
	}
}
