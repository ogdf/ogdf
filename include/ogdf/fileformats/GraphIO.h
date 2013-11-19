/*
 * $Revision: 3831 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 10:00:32 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declares class GraphIO which provides access to all
 *        graph read and write functionality.
 *
 * \author Carsten Gutwenger, Markus Chimani, Karsten Klein, Matthias Woste, Łukasz Hanuszczak
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

#ifndef OGDF_GRAPH_IO_H
#define OGDF_GRAPH_IO_H


#include <ogdf/basic/GridLayout.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraph.h>
#include <sstream>


namespace ogdf {


//! Utility class providing graph I/O in various exchange formats.
class OGDF_EXPORT GraphIO
{
public:
	class SVGSettings
	{
		double m_margin;

		int    m_fontSize;
		string m_fontColor;
		string m_fontFamily;

	public:
		SVGSettings();

		//! Returns the size of the margin around the drawing.
		double margin() const { return m_margin; }

		//! Returns the default font size (font height in pixels).
		int fontSize() const { return m_fontSize; }

		//! Returns the default font color.
		const string &fontColor() const { return m_fontColor; }

		//! Returns the default font family.
		const string &fontFamily() const { return m_fontFamily; }


		//! Sets the size of the margin around the drawing to \a m.
		void margin(double m) { m_margin = m; }

		//! Sets the default font size (font height in pixels) to \a fs.
		void fontSize(int fs) { m_fontSize = fs; }

		//! Sets the default font color to \a fc.
		void fontColor(const string &fc) { m_fontColor = fc; }

		//! Sets the default font family to \a fm.
		void fontFamily(const string &fm) { m_fontFamily = fm; }
	};


	/**
	 * @name Graphs
	 * These functions read and write graphs (instances of type Graph) in various graph formats.
	 */
	//@{

	//! Reads graph \a G in GML format from file \a filename.
	/**
	 * \sa readGML(Graph &G, istream &is) for more details.<br>
	 *     writeGML(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(Graph &G, const char *filename);

	//! Reads graph \a G in GML format from file \a filename.
	/**
	 * \sa readGML(Graph &G, istream &is) for more details.<br>
	 *     writeGML(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(Graph &G, const string &filename);

	//! Reads graph \a G in GML format from input stream \a is.
	/**
	 * The GML (<i>%Graph Modelling Language</i>) file format is an Ascii-based format that has been
	 * developed by Michael Himsolt at the University of Passau. Its full specification can be found in
	 * <a href="http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf">this
	 * technical report</a>.
	 *
	 * \sa writeGML(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(Graph &G, istream &is);

	//! Writes graph \a G in GML format to file \a filename.
	/**
	 * \sa writeGML(const Graph &G, ostream &os) for more details.<br>
	 *     readGML(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const Graph &G, const char *filename);

	//! Writes graph \a G in GML format to file \a filename.
	/**
	 * \sa writeGML(const Graph &G, ostream &os) for more details.<br>
	 *     readGML(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const Graph &G, const string &filename);

	//! Writes graph \a G in GML format to output stream \a os.
	/**
	 * The GML (<i>%Graph Modelling Language</i>) file format is an Ascii-based format that has been
	 * developed by Michael Himsolt at the University of Passau. Its full specification can be found in
	 * <a href="http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf">this
	 * technical report</a>. The GML format stores the basic graph structure, i.e., nodes and edges.
	 *
	 * \sa readGML(Graph &G, istream &is)
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const Graph &G, ostream &os);


	//! Reads graph \a G in OGML format from file \a filename.
	/**
	 * \sa readOGML(Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(Graph &G, const char *filename);

	//! Reads graph \a G in OGML format from file \a filename.
	/**
	 * \sa readOGML(Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(Graph &G, const string &filename);

	//! Reads graph \a G in OGML format from input stream \a is.
	/**
	 * \sa writeOGML(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(Graph &G, istream &is);

	//! Writes graph \a G in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const Graph &G, ostream &os) for more details.<br>
	 *     readOGML(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const Graph &G, const char *filename);

	//! Writes graph \a G in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const Graph &G, ostream &os) for more details.<br>
	 *     readOGML(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const Graph &G, const string &filename);

	//! Writes graph \a G in OGML format to output stream \a os.
	/**
	 * \sa bool readOGML(Graph &G, istream &is)
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const Graph &G, ostream &os);


	//! Reads graph \a G in Rome-Lib format from file \a filename.
	/**
	 * \sa readRome(Graph &G, istream &is) for more details.<br>
	 *     writeRome(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRome(Graph &G, const char *filename);

	//! Reads graph \a G in Rome-Lib format from file \a filename.
	/**
	 * \sa readRome(Graph &G, istream &is) for more details.<br>
	 *     writeRome(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRome(Graph &G, const string &filename);

	//! Reads graph \a G in Rome-Lib format from input stream \a is.
	/**
	 * The Rome-Lib format contains n "node-lines", 1 "separator-line", m "edge-lines" (in this order).
	 * These lines are as follows (whereby all IDs are integer numbers):
	 *  - <b>node-line:</b> <i>NodeId</i> <tt>0</TT>
	 *  - <b>separator-line:</b> starts with a <tt>#</tt>-sign
	 *  - <b>edge-line:</b> <i>EdgeId</i> <tt>0</tt> <i>SourceNodeId</i> <i>TargetNodeId</i>
	 *
	 * \sa writeRome(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRome(Graph &G, istream &is);

	//! Writes graph \a G in Rome-Lib format to file \a filename.
	/**
	 * \sa writeRome(const Graph &G, ostream &os) for more details.<br>
	 *     readRome(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRome(const Graph &G, const char *filename);

	//! Writes graph \a G in Rome-Lib format to file \a filename.
	/**
	 * \sa writeRome(const Graph &G, ostream &os) for more details.<br>
	 *     readRome(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRome(const Graph &G, const string &filename);

	//! Writes graph \a G in Rome-Lib format to output stream \a os.
	/**
	 * \sa readRome(Graph &G, istream &is) for more details about the format.
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRome(const Graph &G, ostream &os);


	//! Reads graph \a G in LEDA graph format from file \a filename.
	/**
	 * \sa readLEDA(Graph &G, istream &is) for more details.<br>
	 *     writeLEDA(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readLEDA(Graph &G, const char *filename);

	//! Reads graph \a G in LEDA graph format from file \a filename.
	/**
	 * \sa readLEDA(Graph &G, istream &is) for more details.<br>
	 *     writeLEDA(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readLEDA(Graph &G, const string &filename);

	//! Reads graph \a G in LEDA graph format from input stream \a is.
	/**
	 * The LEDA graph format is a simple, Ascii-based file format used by the
	 * <a href="http://www.algorithmic-solutions.com/leda/">LEDA library</a>.
	 * Its specification is described in the
	 * <a href="http://www.algorithmic-solutions.info/leda_guide/graphs/leda_native_graph_fileformat.html">
	 * LEDA Guide</a>.
	 *
	 * \sa  writeLEDA(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readLEDA(Graph &G, istream &is);


	//! Writes graph \a G in LEDA graph format to file \a filename.
	/**
	 * \sa writeLEDA(const Graph &G, ostream &os) for more details.<br>
	 *     readLEDA(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeLEDA(const Graph &G, const char *filename);

	//! Writes graph \a G in LEDA graph format to file \a filename.
	/**
	 * \sa writeLEDA(const Graph &G, ostream &os) for more details.<br>
	 *     readLEDA(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeLEDA(const Graph &G, const string &filename);

	//! Writes graph \a G in LEDA graph format to output stream \a os.
	/**
	 * The LEDA graph format is a simple, Ascii-based file format used by the
	 * <a href="http://www.algorithmic-solutions.com/leda/">LEDA library</a>.
	 * Its specification is described in the
	 * <a href="http://www.algorithmic-solutions.info/leda_guide/graphs/leda_native_graph_fileformat.html">
	 * LEDA Guide</a>.
	 *
	 * \sa readLEDA(Graph &G, istream &is)
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeLEDA(const Graph &G, ostream &os);


	//! Reads graph \a G in Chaco format from file \a filename.
	/**
	 * \sa readChaco(Graph &G, istream &is) for more details.<br>
	 *     writeChaco(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readChaco(Graph &G, const char *filename);

	//! Reads graph \a G in Chaco format from file \a filename.
	/**
	 * \sa readChaco(Graph &G, istream &is) for more details.<br>
	 *     writeChaco(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readChaco(Graph &G, const string &filename);

	//! Reads graph \a G in Chaco format from input stream \a is.
	/**
	 * This simple graph format is used by graph partitioning tools like
	 * Chaco, Metis, or Jostle.
	 * Its specification is described in the
	 * <a href="http://staffweb.cms.gre.ac.uk/~wc06/jostle/jostle-exe.pdf">
	 * Jostle User Guide</a>.
	 *
	 * \sa writeChaco(const Graph &G, ostream &os)
	 *
	 * @param G  is assigned the read graph.
	 * @param is is the input stream to be read.
	 * \return true if successful, false otherwise.
	 * */
	static bool readChaco(Graph &G, istream &is);

	//! Writes graph \a G in Chaco format to file \a filename.
	/**
	 * \sa writeChaco(const Graph &G, ostream &os) for more details.<br>
	 *     readChaco(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeChaco(const Graph &G, const char *filename);

	//! Writes graph \a G in Chaco format to file \a filename.
	/**
	 * \sa writeChaco(const Graph &G, ostream &os) for more details.<br>
	 *     readChaco(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeChaco(const Graph &G, const string &filename);

	//! Writes graph \a G in Chaco format to output stream \a os.
	/**
	 * This simple graph format is used by graph partitioning tools like
	 * Chaco, Metis, or Jostle.
	 * Its specification is described in the
	 * <a href="http://staffweb.cms.gre.ac.uk/~wc06/jostle/jostle-exe.pdf">
	 * Jostle User Guide</a>.
	 *
	 * \sa readChaco(Graph &G, istream &is)
	 *
	 * @param G  is the graph to be written.
	 * @param os is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeChaco(const Graph &G, ostream &os);


	//! Reads graph \a G in a simple format as used in Petra Mutzel's thesis from file \a filename.
	/**
	 * \sa readPMDissGraph(Graph &G, istream &is) for more details.<br>
	 *     writePMDissGraph(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readPMDissGraph(Graph &G, const char *filename);

	//! Reads graph \a G in a simple format as used in Petra Mutzel's thesis from file \a filename.
	/**
	 * \sa readPMDissGraph(Graph &G, istream &is) for more details.<br>
	 *     writePMDissGraph(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readPMDissGraph(Graph &G, const string &filename);

	//! Reads graph \a G in a simple format as used in Petra Mutzel's thesis from input stream \a is.
	/**
	 * This simple graph format has a leading line stating the name of the graph
	 * and a following line stating the size of the graph:
	 *
	 * <pre>
	 * *BEGIN unknown_name.numN.numE
	 * *GRAPH numN numE UNDIRECTED UNWEIGHTED
	 * </pre>
	 *
	 * \sa writePMDissGraph(const Graph &G, ostream &os)
	 *
	 * @param G  is assigned the read graph.
	 * @param is is the input stream to be read.
	 * \return true if successful, false otherwise.
	 * */
	static bool readPMDissGraph(Graph &G, istream &is);

	//! Writes graph \a G in a simple format as used in Petra Mutzel's thesis to file \a filename.
	/**
	 * \sa writePMDissGraph(const Graph &G, ostream &os) for more details.<br>
	 *     readPMDissGraph(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writePMDissGraph(const Graph &G, const char *filename);

	//! Writes graph \a G in a simple format as used in Petra Mutzel's thesis to file \a filename.
	/**
	 * \sa writePMDissGraph(const Graph &G, ostream &os) for more details.<br>
	 *     readPMDissGraph(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writePMDissGraph(const Graph &G, const string &filename);

	//! Writes graph \a G in a simple format as used in Petra Mutzel's thesis to output stream \a os.
	/**
	 * This simple graph format has a leading line stating the name of the graph
	 * and a following line stating the size of the graph:
	 *
	 * <pre>
	 * *BEGIN unknown_name.numN.numE
	 * *GRAPH numN numE UNDIRECTED UNWEIGHTED
	 * </pre>
	 *
	 * \sa readPMDissGraph(Graph &G, istream &is)
	 *
	 * @param G  is the graph to be written.
	 * @param os is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writePMDissGraph(const Graph &G, ostream &os);


	//! Reads graph \a G in Y-graph format from file \a filename.
	/**
	 * \sa readYGraph(Graph &G, istream &is) for more details.
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readYGraph(Graph &G, const char *filename);

	//! Reads graph \a G in Y-graph format from file \a filename.
	/**
	 * \sa readYGraph(Graph &G, istream &is) for more details.
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readYGraph(Graph &G, const string &filename);

	//! Reads graph \a G in Y-graph format from input stream \a is.
	/**
	 * This format is e.g. produced by NAUTY (http://www.cs.sunysb.edu/~algorith/implement/nauty/implement.shtml).
	 *
	 * Details  on the format, as given in NAUTYs graph generator (see above link):
	 * "[A] graph occupies one line with a terminating newline.
	 * Except for the newline, each byte has the format  01xxxxxx, where
	 * each "x" represents one bit of data.
	 *
	 * First byte:  xxxxxx is the number of vertices n
	 *
	 * Other ceiling(n(n-1)/12) bytes:  These contain the upper triangle of
	 * the adjacency matrix in column major order.  That is, the entries
	 * appear in the order (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(0,4),... .
	 * The bits are used in left to right order within each byte.
	 * Any unused bits on the end are set to zero.
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readYGraph(Graph &G, istream &is);


	//@}
	/**
	 * @name Clustered graphs
	 * These functions read and write clustered graphs (instances of type ClusterGraph)
	 * in various graph formats. Read functions take a pair (\a C, \a G) as parameters,
	 * where \a G is the graph associated with the clustered graph \a C.
	 */
	//@{

	//! Reads clustered graph (\a C, \a G) in GML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in GML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in GML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeGML(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraph &C, Graph &G, istream &is);

	//! Writes clustered graph \a C in GML format to file \a filename.
	/**
	 * \sa writeGML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGML(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in GML format to file \a filename.
	/**
	 * \sa writeGML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGML(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in GML format to output stream \a os.
	/**
	 * \sa readGML(ClusterGraph &C, Graph &G, istream &is)
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraph &C, ostream &os);


	//! Reads clustered graph (\a C, \a G) in OGML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readOGML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in OGML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readOGML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in OGML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeOGML(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraph &C, Graph &G, istream &is);


	//! Writes clustered graph \a C in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readOGML(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readOGML(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in OGML format to output stream \a os.
	/**
	 * \sa readOGML(ClusterGraph &C, Graph &G, istream &is)
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraph &C, ostream &os);


	//@}
	/**
	 * @name Graphs with attributes
	 * These functions read and write graphs with additional attributes (instances of type
	 * GraphAttributes) in various graph formats. Read functions take a pair (\a A, \a G) as parameters,
	 * where \a G is the graph associated with the graph attributes \a A.
	 */
	//@{

	//! Reads graph \a G with attributes \a A in GML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in GML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in GML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeGML(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(GraphAttributes &A, Graph &G, istream &is);

	//! Writes graph with attributes \a A in GML format to file \a filename.
	/**
	 * \sa writeGML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGML(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GML format to file \a filename.
	/**
	 * \sa writeGML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGML(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GML format to output stream \a os.
	/**
	 * \sa readGML(GraphAttributes &A, Graph &G, istream &is)
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const GraphAttributes &A, ostream &os);


	//! Reads graph \a G with attributes \a A in OGML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readOGML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in OGML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readOGML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in OGML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeOGML(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(GraphAttributes &A, Graph &G, istream &is);

	//! Writes graph with attributes \a A in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readOGML(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readOGML(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in OGML format to output stream \a os.
	/**
	 * \sa readOGML(GraphAttributes &A, Graph &G, istream &is)
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const GraphAttributes &A, ostream &os);


	//! Reads graph \a G with edge weights stored in \a A in Rudy format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readRudy(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeRudy(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes (only edge weights (as doubles) are used).
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRudy(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with edge weights stored in \a A in Rudy format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readRudy(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeRudy(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes (only edge weights (as doubles) are used).
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRudy(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with edge weights stored in \a A in Rudy format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeRudy(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes (only edge weights (as doubles) are used).
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readRudy(GraphAttributes &A, Graph &G, istream &is);

	//! Writes graph with edge weights stored in \a A in Rudy format to file \a filename.
	/**
	 * \sa writeRudy(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readRudy(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written (only edge weights
	 *                 (as doubles) are stored in this format).
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRudy(const GraphAttributes &A, const char *filename);

	//! Writes graph with edge weights stored in \a A in Rudy format to file \a filename.
	/**
	 * \sa writeRudy(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readRudy(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written (only edge weights
	 *                 (as doubles) are stored in this format).
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRudy(const GraphAttributes &A, const string &filename);

	//! Writes graph with edge weights stored in \a A in Rudy format to output stream \a os.
	/**
	 * \sa readRudy(GraphAttributes &A, Graph &G, istream &is)
	 *
	 * @param A   specifies the graph and its attributes to be written (only edge weights
	 *            (as doubles) are stored in this format).
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeRudy(const GraphAttributes &A, ostream &os);


	//@}
	/**
	 * @name Clustered graphs with attributes
	 * These functions read and write clustered graphs with additional attributes (instances of type
	 * ClusterGraphAttributes) in various graph formats. Read functions take a triple (\a A, \a C, \a G)
	 * as parameters, where \a C is the clustered graph associated with the graph attributes \a A and
	 * \a G is the graph associated with \a C.
	 */
	//@{

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGML(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GML format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeGML(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);

	//! Writes with attributes \a A in GML format to file \a filename.
	/**
	 * \sa writeGML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GML format to file \a filename.
	/**
	 * \sa writeGML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GML format to output stream \a os.
	/**
	 * \sa readGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGML(const ClusterGraphAttributes &A, ostream &os);


	//! Reads clustered graph (\a C, \a G) with attributes \a A in OGML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in OGML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeOGML(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in OGML format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeOGML(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);

	//! Writes with attributes \a A in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in OGML format to file \a filename.
	/**
	 * \sa writeOGML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in OGML format to output stream \a os.
	/**
	 * \sa readOGML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is)
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeOGML(const ClusterGraphAttributes &A, ostream &os);


	//@}
	/**
	 * @name Hypergraphs as point-based expansion
	 * These functions load hypergraphs stored in file formats used for electrical circuits.
	 * The hypergraphs are directly transformed into their point-based expansions (and hence stored
	 * in a usual Graph and not a Hypergraph).
	 */
	//@{

	//!  Reads a hypergraph (as point-based expansion) in BENCH format from file \a filename.
	/**
	 * A hypergraph in OGDF is represented by its point-based expansion, i.e., for each
	 * hyperedge <i>h</i> we have a corresponding hypernode <i>n</i>. All nodes originally
	 * incident to <i>h</i> are incident to <i>n</i>, i.e., have regular edges to <i>n</i>.
	 *
	 * \warning
	 * This is a very simple implementation only usable for very properly formatted files!
	 *
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the BENCH-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param filename   is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readBENCH(Graph &G, List<node>& hypernodes, List<edge> *shell, const char *filename);

	//!  Reads a hypergraph (as point-based expansion) in BENCH format from file \a filename.
	/**
	 * A hypergraph in OGDF is represented by its point-based expansion, i.e., for each
	 * hyperedge <i>h</i> we have a corresponding hypernode <i>n</i>. All nodes originally
	 * incident to <i>h</i> are incident to <i>n</i>, i.e., have regular edges to <i>n</i>.
	 *
	 * \warning
	 * This is a very simple implementation only usable for very properly formatted files!
	 *
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the BENCH-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param filename   is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readBENCH(Graph &G, List<node>& hypernodes, List<edge> *shell, const string &filename);

	//!  Reads a hypergraph (as point-based expansion) in BENCH format from input stream \a is.
	/**
	 * A hypergraph in OGDF is represented by its point-based expansion, i.e., for each
	 * hyperedge <i>h</i> we have a corresponding hypernode <i>n</i>. All nodes originally
	 * incident to <i>h</i> are incident to <i>n</i>, i.e., have regular edges to <i>n</i>.
	 *
	 * \warning
	 * This is a very simple implementation only usable for very properly formatted files!
	 *
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the BENCH-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param is         is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readBENCH(Graph &G, List<node>& hypernodes, List<edge> *shell, istream &is);


	//! Reads a hypergraph (as point-based expansion) in PLA format from file \a filename.
	/**
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the PLA-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param filename   is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 *
	 * \sa readPLA(Graph &G, List<node>& hypernodes, List<edge> *shell, istream &is) for details
	 */
	static bool readPLA(Graph &G, List<node>& hypernodes, List<edge>* shell, const char *filename);

	//! Reads a hypergraph (as point-based expansion) in PLA format from file \a filename.
	/**
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the PLA-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param filename   is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 *
	 * \sa readPLA(Graph &G, List<node>& hypernodes, List<edge> *shell, istream &is) for details
	 */
	static bool readPLA(Graph &G, List<node>& hypernodes, List<edge>* shell, const string &filename);

	//!  Reads a hypergraph (as point-based expansion) in PLA format from input stream \a is.
	/**
	 * A hypergraph in OGDF is represented by its point-based expansion, i.e., for each
	 * hyperedge <i>h</i> we have a corresponding hypernode <i>n</i>. All nodes originally
	 * incident to <i>h</i> are incident to <i>n</i>, i.e., have regular edges to <i>n</i>.
	 *
	 * \warning
	 * This is a very simple implementation only usable for very properly formatted files!
	 *
	 * @param G          is assigned the read graph (point-based expansion of the hypergraph).
	 * @param hypernodes is assigned the list of nodes which have to be interpreted as hypernodes.
	 * @param shell      if 0 only the PLA-hypergraph is read. Otherwise we extend the graph
	 *                   by a simple edge <i>e=(i,o)</i> and two hyperedges: one hyperedges groups all input nodes and
	 *                   <i>i</i> together, the other hyperedge groups all output edges and <i>o</i>.
	 *                   These additional edges are then also collocated in shell.
	 * @param is         is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readPLA(Graph &G, List<node>& hypernodes, List<edge> *shell, istream &is);


	//@}
	/**
	 * @name Graphs Drawing Contest formats
	 * These functions read and write graphs in a simple text-based file formats that have been used
	 * in the Graph Drawing Contest.
	 */
	//@{

	//! Reads graph \a G with grid layout \a gl in GD-Challenge-format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa writeChallengeGraph(const Graph &G, const GridLayout &gl, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param gl       is assigned the grid layout.
	 * @param filename is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readChallengeGraph(Graph &G, GridLayout &gl, const char *filename);

	//! Reads graph \a G with grid layout \a gl in GD-Challenge-format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa writeChallengeGraph(const Graph &G, const GridLayout &gl, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param gl       is assigned the grid layout.
	 * @param filename is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readChallengeGraph(Graph &G, GridLayout &gl, const string &filename);

	//! Reads graph \a G with grid layout \a gl in GD-Challenge-format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa writeChallengeGraph(const Graph &G, const GridLayout &gl, ostream &os)
	 *
	 * @param G  is assigned the read graph.
	 * @param gl is assigned the grid layout.
	 * @param is is the input stream from which the graph is read.
	 * \return true if successful, false otherwise.
	 */
	static bool readChallengeGraph(Graph &G, GridLayout &gl, istream &is);

	//! Writes graph \a G with grid layout \a gl in GD-Challenge-format to output stream \a os.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa readChallengeGraph(Graph &G, GridLayout &gl, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param gl       specifies the grid layout of \a G to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeChallengeGraph(const Graph &G, const GridLayout &gl, const char *filename);

	//! Writes graph \a G with grid layout \a gl in GD-Challenge-format to output stream \a os.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa readChallengeGraph(Graph &G, GridLayout &gl, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param gl       specifies the grid layout of \a G to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeChallengeGraph(const Graph &G, const GridLayout &gl, const string &filename);

	//! Writes graph \a G with grid layout \a gl in GD-Challenge-format to output stream \a os.
	/**
	 * \pre \a G is the graph associated with grid layout \a gl.
	 * \sa readChallengeGraph(Graph &G, GridLayout &gl, istream &is)
	 *
	 * @param G  is the graph to be written.
	 * @param gl specifies the grid layout of \a G to be written.
	 * @param os is the output stream to which the graph is written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeChallengeGraph(const Graph &G, const GridLayout &gl, ostream &os);

	//! Reads graph \a G in GraphML format from file \a filename.
	/**
	 * \sa readGraphML(Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(Graph &G, const char *filename);

	//! Reads graph \a G in GraphML format from file \a filename.
	/**
	 * \sa readGraphML(Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(Graph &G, const string &filename);

	//! Reads graph \a G in GraphML format from input stream \a is.
	/**
	 * \sa writeGraphML(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) in GraphML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGraphML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in GraphML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGraphML(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in GraphML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeGraphML(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraph &C, Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in GraphML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGraphML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in GraphML format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGraphML(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in GraphML format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeGraphML(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(GraphAttributes &A, Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GraphML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GraphML format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGraphML(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GraphML format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeGraphML(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);

	//! Writes graph \a G in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const Graph &G, ostream &os) for more details.<br>
	 *     readGraphML(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const Graph &G, const char *filename);

	//! Writes graph \a G in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const Graph &G, ostream &os) for more details.<br>
	 *     readGraphML(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const Graph &G, const string &filename);

	//! Writes graph \a G in GraphML format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const Graph &G, ostream &os);

	//! Writes clustered graph \a C in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGraphML(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGraphML(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in GraphML format to output stream \a os.
	/**
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraph &C, ostream &os);

	//! Writes graph with attributes \a A in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGraphML(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGraphML(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GraphML format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const GraphAttributes &A, ostream &os);

	//! Writes with attributes \a A in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GraphML format to file \a filename.
	/**
	 * \sa writeGraphML(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGraphML(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GraphML format to output stream \a os.
	/**
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGraphML(const ClusterGraphAttributes &A, ostream &os);

	//! Reads graph \a G in DOT format from file \a filename.
	/**
	 * \sa readDOT(Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(Graph &G, const char *filename);

	//! Reads graph \a G in DOT format from file \a filename.
	/**
	 * \sa readDOT(Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(Graph &G, const string &filename);

	//! Reads graph \a G in DOT format from input stream \a is.
	/**
	 * \sa writeDOT(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) in DOT format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readDOT(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in DOT format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readDOT(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in DOT format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeDOT(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraph &C, Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in DOT format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readDOT(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in DOT format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readDOT(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in DOT format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeDOT(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(GraphAttributes &A, Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in DOT format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in DOT format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeDOT(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in DOT format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeDOT(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);

	//! Writes graph \a G in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const Graph &G, ostream &os) for more details.<br>
	 *     readDOT(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const Graph &G, const char *filename);

	//! Writes graph \a G in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const Graph &G, ostream &os) for more details.<br>
	 *     readDOT(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const Graph &G, const string &filename);

	//! Writes graph \a G in DOT format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const Graph &G, ostream &os);

	//! Writes clustered graph \a C in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readDOT(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readDOT(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in DOT format to output stream \a os.
	/**
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraph &C, ostream &os);

	//! Writes graph with attributes \a A in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readDOT(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readDOT(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in DOT format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const GraphAttributes &A, ostream &os);

	//! Writes with attributes \a A in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in DOT format to file \a filename.
	/**
	 * \sa writeDOT(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readDOT(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in DOT format to output stream \a os.
	/**
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDOT(const ClusterGraphAttributes &A, ostream &os);

	//! Reads graph \a G in GEXF format from file \a filename.
	/**
	 * \sa readGEXF(Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(Graph &G, const char *filename);

	//! Reads graph \a G in GEXF format from file \a filename.
	/**
	 * \sa readGEXF(Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(Graph &G, const string &filename);

	//! Reads graph \a G in GEXF format from input stream \a is.
	/**
	 * \sa writeGEXF(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) in GEXF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGEXF(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in GEXF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readGEXF(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in GEXF format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeGEXF(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraph &C, Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in GEXF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGEXF(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in GEXF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGEXF(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in GEXF format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeGEXF(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(GraphAttributes &A, Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GEXF format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GEXF format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeGEXF(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in GEXF format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeGEXF(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);

	//! Writes graph \a G in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const Graph &G, ostream &os) for more details.<br>
	 *     readGEXF(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const Graph &G, const char *filename);

	//! Writes graph \a G in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const Graph &G, ostream &os) for more details.<br>
	 *     readGEXF(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const Graph &G, const string &filename);

	//! Writes graph \a G in GEXF format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const Graph &G, ostream &os);

	//! Writes clustered graph \a C in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGEXF(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readGEXF(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in GEXF format to output stream \a os.
	/**
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraph &C, ostream &os);

	//! Writes graph with attributes \a A in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGEXF(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGEXF(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GEXF format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const GraphAttributes &A, ostream &os);

	//! Writes with attributes \a A in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GEXF format to file \a filename.
	/**
	 * \sa writeGEXF(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readGEXF(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GEXF format to output stream \a os.
	/**
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGEXF(const ClusterGraphAttributes &A, ostream &os);

	//! Reads graph \a G in GDF format from file \a filename.
	/**
	 * \sa readGDF(Graph &G, istream &is) for more details.<br>
	 *     writeGDF(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(Graph &G, const char *filename);

	//! Reads graph \a G in GDF format from file \a filename.
	/**
	 * \sa readGDF(Graph &G, istream &is) for more details.<br>
	 *     writeGDF(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(Graph &G, const string &filename);

	//! Reads graph \a G in GDF format from input stream \a is.
	/**
	 * \sa writeGDF(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in GDF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGDF(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGDF(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in GDF format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readGDF(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeGDF(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in GDF format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeGDF(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readGDF(GraphAttributes &A, Graph &G, istream &is);

	//! Writes graph \a G in GDF format to file \a filename.
	/**
	 * \sa writeGDF(const Graph &G, ostream &os) for more details.<br>
	 *     readGDF(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const Graph &G, const char *filename);

	//! Writes graph \a G in GDF format to file \a filename.
	/**
	 * \sa writeGDF(const Graph &G, ostream &os) for more details.<br>
	 *     readGDF(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const Graph &G, const string &filename);

	//! Writes graph \a G in GDF format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const Graph &G, ostream &os);

	//! Writes graph with attributes \a A in GDF format to file \a filename.
	/**
	 * \sa writeGDF(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGDF(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in GDF format to file \a filename.
	/**
	 * \sa writeGDF(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readGDF(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in GDF format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeGDF(const GraphAttributes &A, ostream &os);

	//! Reads graph \a G in TLP format from file \a filename.
	/**
	 * \sa readTLP(Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(Graph &G, const char *filename);

	//! Reads graph \a G in TLP format from file \a filename.
	/**
	 * \sa readTLP(Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(Graph &G, const string &filename);

	//! Reads graph \a G in TLP format from input stream \a is.
	/**
	 * \sa writeTLP(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) in TLP format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readTLP(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const ClusterGraph &C, const char *filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) in TLP format from file \a filename.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa readTLP(ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const ClusterGraph &C, const string &filename)
	 *
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) in TLP format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with clustered graph \a C.
	 * \sa writeTLP(const ClusterGraph &C, ostream &os)
	 *
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraph &C, Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in TLP format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readTLP(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in TLP format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readTLP(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in TLP format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeTLP(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(GraphAttributes &A, Graph &G, istream &is);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in TLP format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const ClusterGraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in TLP format from file \a filename.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is) for more details.<br>
	 *     writeTLP(const ClusterGraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param C        is assigned the read clustered graph (cluster structure).
	 * @param G        is assigned the read clustered graph (graph structure).
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename);

	//! Reads clustered graph (\a C, \a G) with attributes \a A in TLP format from input stream \a is.
	/**
	 * \pre \a C is the clustered graph associated with attributes \a A, and \a G is the graph associated with \a C.
	 * \sa writeTLP(const ClusterGraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param C   is assigned the read clustered graph (cluster structure).
	 * @param G   is assigned the read clustered graph (graph structure).
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, istream &is);


	//! Writes graph \a G in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const Graph &G, ostream &os) for more details.<br>
	 *     readTLP(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const Graph &G, const char *filename);

	//! Writes graph \a G in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const Graph &G, ostream &os) for more details.<br>
	 *     readTLP(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const Graph &G, const string &filename);

	//! Writes graph \a G in TLP format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const Graph &G, ostream &os);

	//! Writes clustered graph \a C in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readTLP(ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraph &C, const char *filename);

	//! Writes clustered graph \a C in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const ClusterGraph &C, ostream &os) for more details.<br>
	 *     readTLP(ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param C        is the clustered graph to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraph &C, const string &filename);

	//! Writes clustered graph \a C in TLP format to output stream \a os.
	/**
	 *
	 * @param C   is the clustered graph to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraph &C, ostream &os);

	//! Writes graph with attributes \a A in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readTLP(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readTLP(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in TLP format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const GraphAttributes &A, ostream &os);

	//! Writes with attributes \a A in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const char *filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in TLP format to file \a filename.
	/**
	 * \sa writeTLP(const ClusterGraphAttributes &A, ostream &os) for more details.<br>
	 *     readTLP(ClusterGraphAttributes &A, ClusterGraph &C, Graph &G, const string &filename)
	 *
	 * @param A        specifies the clustered graph and its attributes to be written.
	 * @param filename is the name of the file to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in TLP format to output stream \a os.
	/**
	 *
	 * @param A   specifies the clustered graph and its attributes to be written.
	 * @param os  is the output stream to which the clustered graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeTLP(const ClusterGraphAttributes &A, ostream &os);

	//! Reads graph \a G in DL format from file \a filename.
	/**
	 * \sa readDL(Graph &G, istream &is) for more details.<br>
	 *     writeDL(const Graph &G, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(Graph &G, const char *filename);

	//! Reads graph \a G in DL format from file \a filename.
	/**
	 * \sa readDL(Graph &G, istream &is) for more details.<br>
	 *     writeDL(const Graph &G, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(Graph &G, const string &filename);

	//! Reads graph \a G in DL format from input stream \a is.
	/**
	 * \sa writeDL(const Graph &G, ostream &os)
	 *
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(Graph &G, istream &is);

	//! Reads graph \a G with attributes \a A in DL format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readDL(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeDL(const GraphAttributes &A, const char *filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(GraphAttributes &A, Graph &G, const char *filename);

	//! Reads graph \a G with attributes \a A in DL format from file \a filename.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa readDL(GraphAttributes &A, Graph &G, istream &is) for more details.<br>
	 *     writeDL(const GraphAttributes &A, const string &filename)
	 *
	 * @param A        is assigned the graph's attributes.
	 * @param G        is assigned the read graph.
	 * @param filename is the name of the file to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(GraphAttributes &A, Graph &G, const string &filename);

	//! Reads graph \a G with attributes \a A in DL format from input stream \a is.
	/**
	 * \pre \a G is the graph associated with attributes \a A.
	 * \sa writeDL(const GraphAttributes &A, ostream &os)
	 *
	 * @param A   is assigned the graph's attributes.
	 * @param G   is assigned the read graph.
	 * @param is  is the input stream to be read.
	 * @return true if successful, false otherwise.
	 */
	static bool readDL(GraphAttributes &A, Graph &G, istream &is);

	//! Writes graph \a G in DL format to file \a filename.
	/**
	 * \sa writeDL(const Graph &G, ostream &os) for more details.<br>
	 *     readDL(Graph &G, const char *filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const Graph &G, const char *filename);

	//! Writes graph \a G in DL format to file \a filename.
	/**
	 * \sa writeDL(const Graph &G, ostream &os) for more details.<br>
	 *     readDL(Graph &G, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const Graph &G, const string &filename);

	//! Writes graph \a G in DL format to output stream \a os.
	/**
	 *
	 * @param G   is the graph to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const Graph &G, ostream &os);

	//! Writes graph with attributes \a A in DL format to file \a filename.
	/**
	 * \sa writeDL(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readDL(GraphAttributes &A, Graph &G, const char *filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const GraphAttributes &A, const char *filename);

	//! Writes graph with attributes \a A in DL format to file \a filename.
	/**
	 * \sa writeDL(const GraphAttributes &A, ostream &os) for more details.<br>
	 *     readDL(GraphAttributes &A, Graph &G, const string &filename)
	 *
	 * @param A        specifies the graph and its attributes to be written.
	 * @param filename is the name of the file to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const GraphAttributes &A, const string &filename);

	//! Writes graph with attributes \a A in DL format to output stream \a os.
	/**
	 *
	 * @param A   specifies the graph and its attributes to be written.
	 * @param os  is the output stream to which the graph will be written.
	 * @return true if successful, false otherwise.
	 */
	static bool writeDL(const GraphAttributes &A, ostream &os);

	//@}
	/**
	 * @name SteinLib instances
	 * These functions read SteinLib instances stored in STP format and convert them into a weighted graph (represented by
	 * an EdgeWeightedGraph) and a list of terminal nodes.
	 */
	//@{

	//! Reads a SteinLib instance from file \a filename and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \sa readSTP(EdgeWeightedGraph<double> &wG, List<node> &terminals, NodeArray<bool> &isTerminal, istream &is) for more details.
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  filename   is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<double> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		const char           *filename);

	//! Reads a SteinLib instance from file \a filename and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \sa readSTP(EdgeWeightedGraph<int> &wG, List<node> &terminals, NodeArray<bool> &isTerminal, istream &is) for more details.
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  filename   is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<int> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		const char           *filename);

	//! Reads a SteinLib instance from file \a filename and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \sa readSTP(EdgeWeightedGraph<double> &wG, List<node> &terminals, NodeArray<bool> &isTerminal, istream &is) for more details.
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  filename   is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<double> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		const string         &filename);

	//! Reads a SteinLib instance from file \a filename and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \sa readSTP(EdgeWeightedGraph<int> &wG, List<node> &terminals, NodeArray<bool> &isTerminal, istream &is) for more details.
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  filename   is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<int> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		const string         &filename);

	//! Reads a SteinLib instance from input stream \a is and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \warning The coordinate section of the SteinLib instance is not read!
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  is is the input stream from which the graph is read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<double> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		istream              &is);

	//! Reads a SteinLib instance from input stream \a is and converts it into a weighted graph \a wG and a set of terminal nodes \a terminals.
	/**
	 * \warning The coordinate section of the SteinLib instance is not read!
	 *
	 * @param  wG         is assigned the graph stored in the SteinLib instance.
	 * @param  terminals  is assgined the list of terminals specified in the SteinLib instance.
	 * @param  isTerminal is assigned the incidence vector for the terminal nodes.
	 * @param  is is the input stream from which the graph is read.
	 * \return true if successful, false otherwise.
	 */
	static bool readSTP(
		EdgeWeightedGraph<int> &wG,
		List<node>           &terminals,
		NodeArray<bool>      &isTerminal,
		istream              &is);


	//@}
	/**
	 * @name Graphs with subgraph
	 * These functions read and write graphs in a simple text-based file format that also specifies
	 * a subgraph (given as a list of edges).
	 */
	//@{

	//! Reads graph \a G with subgraph defined by \a delEdges from file \a filename.
	/**
	 * \sa writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const char *filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param delEdges is assigned the edges of the subgraph.
	 * @param filename is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const char *filename);

	//! Reads graph \a G with subgraph defined by \a delEdges from file \a filename.
	/**
	 * \sa writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const string &filename)
	 *
	 * @param G        is assigned the read graph.
	 * @param delEdges is assigned the edges of the subgraph.
	 * @param filename is the name of the file to be read.
	 * \return true if successful, false otherwise.
	 */
	static bool readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const string &filename);

	//! Reads graph \a G with subgraph defined by \a delEdges from stream \a is.
	/**
	 * \sa writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, ostream &os)
	 *
	 * @param G        is assigned the read graph.
	 * @param delEdges is assigned the edges of the subgraph.
	 * @param          is is the input stream from which the graph is read.
	 * \return true if successful, false otherwise.
	 */
	static bool readEdgeListSubgraph(Graph &G, List<edge> &delEdges, istream &is);

	//! Writes graph \a G with subgraph defined by \a delEdges to file \a filename.
	/**
	 * \sa readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const char *filename)
	 *
	 * @param G is the graph to be written.
	 * @param delEdges specifies the edges of the subgraph to be stored.
	 * @param filename is the name of the file to which the graph will be written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const char *filename);

	//! Writes graph \a G with subgraph defined by \a delEdges to file \a filename.
	/**
	 * \sa readEdgeListSubgraph(Graph &G, List<edge> &delEdges, const string &filename)
	 *
	 * @param G        is the graph to be written.
	 * @param delEdges specifies the edges of the subgraph to be stored.
	 * @param filename is the name of the file to which the graph will be written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, const string &filename);

	//! Writes graph \a G with subgraph defined by \a delEdges to stream \a os.
	/**
	 * \sa readEdgeListSubgraph(Graph &G, List<edge> &delEdges, istream &is)
	 *
	 * @param G        is the graph to be written.
	 * @param delEdges specifies the edges of the subgraph to be stored.
	 * @param os       is the output stream to which the graph will be written.
	 * \return true if successful, false otherwise.
	 */
	static bool writeEdgeListSubgraph(const Graph &G, const List<edge> &delEdges, ostream &os);


	//@}
	/**
	 * @name Graphics formats
	 * These functions draw graphs and export them as SVG (Scalable Vector Graphics) vectors graphics.
	 */
	//@{

	static bool drawSVG(const GraphAttributes &A, const char *filename, const SVGSettings &settings = svgSettings);
	static bool drawSVG(const GraphAttributes &A, const string &filename, const SVGSettings &settings = svgSettings);
	static bool drawSVG(const GraphAttributes &A, ostream &os, const SVGSettings &settings = svgSettings);

	static bool drawSVG(const ClusterGraphAttributes &A, const char *filename, const SVGSettings &settings = svgSettings);
	static bool drawSVG(const ClusterGraphAttributes &A, const string &filename, const SVGSettings &settings = svgSettings);
	static bool drawSVG(const ClusterGraphAttributes &A, ostream &os, const SVGSettings &settings = svgSettings);


	//@}
	/**
	 * @name Utility functions for indentation
	 * Text based write methods that use indentation for better readability of the produced text files
	 * apply a customizable indentation character and indentation width.
	 */
	//@{

	//! Returns the currently used indentation character.
	static char indentChar() { return s_indentChar; }

	//! Returns the currently used indentation width.
	static int indentWidth() { return s_indentWidth; }

	//! Sets the indentation character to \ c.
	/**
	 * \pre \a c must be a white-space character (e.g., a space or a tab).
	 */
	static void setIndentChar(char c) {
		OGDF_ASSERT(isspace((int)c));
		s_indentChar = c;
	}

	//! Sets the indentation width to \a w.
	/**
	 * \pre \a w must be non-negative.
	 * Setting the indentation width to 0 suppresses indentation.
	 */
	static void setIndentWidth(int w) {
		if(w >= 0) s_indentWidth = w;
	}

	//! Prints indentation for indentation \a depth to output stream \a os and returns \a os.
	static ostream &indent(ostream &os, int depth);

	//@}


	static SVGSettings svgSettings;

private:
	static char s_indentChar;	//!< Character used for indentation.
	static int  s_indentWidth;	//!< Number of indent characters used for indentation.
};


} // end namespace ogdf


#endif
