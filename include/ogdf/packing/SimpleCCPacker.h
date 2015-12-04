/** \file
* \brief Splits and packs the components of a Graph.
*
* \brief Simple proxy class that uses the TileToRowsCCPacker.
*        Use it for layouts that do not support disconnected graphs.
*
* \author Martin Gronemann
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

#pragma once

#include <ogdf/module/LayoutModule.h>


namespace ogdf {

//! A Layout Module for using the TileToRowsCCPacker
class SimpleCCPackerModule : public LayoutModule
{
public:
    // !Constructor for the CCPacker
    SimpleCCPackerModule(LayoutModule* pSubLayoutModule = nullptr) : m_pSubLayoutModule(pSubLayoutModule)
    {
        m_leftMargin = 10.0;
        m_rightMargin = 10.0;
        m_bottomMargin = 10.0;
        m_topMargin = 10.0;
    }

    void setMargins(double left, double top, double right, double bottom)
    {
        m_leftMargin = left;
        m_rightMargin = right;
        m_bottomMargin = bottom;
        m_topMargin = top;
    }

    virtual void call(GraphAttributes& GA) override;

protected:
    double m_leftMargin;
    double m_rightMargin;
    double m_bottomMargin;
    double m_topMargin;

    void computeBoundingBox(const GraphAttributes& graphAttributes, DPoint& min_coord, DPoint& max_coord );
    LayoutModule* m_pSubLayoutModule;
};

} // end of namespace ogdf
