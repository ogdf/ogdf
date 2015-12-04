/** \file
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

#pragma once

#include <ogdf/basic/basic.h>

namespace ogdf {

template<int Dim>
inline double computeDeltaAndDistance(const double a[Dim], const double b[Dim], double delta[Dim])
{
    // distance var
    double dist = 0;

    // for all dim
    for (int d = 0; d < Dim; d++)
    {
        // delta in d dim
        delta[d] = a[d] - b[d];

        // squared distance sum
        dist += delta[d] * delta[d];
    }

    // distance square root
    return sqrt(dist);
}

template<>
inline double computeDeltaAndDistance<2>(const double a[2], const double b[2], double delta[2])
{
    // delta in d dim
    delta[0] = a[0] - b[0];
    delta[1] = a[1] - b[1];

    // distance square root
    return sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
}

template<int Dim, int K>
inline void RepForceFunctionNewton(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;

    double t = dist;
    for (int i = 1; i < K; i++)
    {
        t *= dist;
    }

    force        = 1.0 / (t);
    force_prime  =  (double)(K) / (t * dist);
}

template<>
inline void RepForceFunctionNewton<2,2>(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;
    force_prime  =  2.0 / (dist * dist * dist);
    force        = force_prime * dist * 0.5;
}

template<>
inline void RepForceFunctionNewton<2,1>(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;
    force_prime  =  1.0 / (dist * dist);
    force        =  1.0 / (dist);
}


template<int Dim>
inline void AttrForceFunctionLog(double dist, double &force, double& force_prime)
{
    // distance square root
    force        = log(dist / 1.0);
    force_prime  =  1.0 / dist;
}


template<int Dim, int K>
inline void AttrForceFunctionPow(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;

    double t = dist;
    for (int i = 1; i < K; i++)
    {
        t *= dist;
    }

    force        =  t;
    force_prime  =  (double)(K) * (t / dist);
}

template<>
inline void AttrForceFunctionPow<2,2>(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;

    force        =  dist * dist;
    force_prime  =  2.0 * dist;
}


template<>
inline void AttrForceFunctionPow<2,1>(double dist, double &force, double& force_prime)
{
    // distance square root
    dist = dist + 0.1;

    force        =  dist;
    force_prime  =  1.0;
}

/*

template<int Dim>
inline void RepForceFunctionInvGauss(const double a[Dim], const double b[Dim], double result[Dim])
{
    // the range of the repulsive force
    const double force_range = 10.0;

    // the amount of repulsive force
    const double force_amount = 20.0;

    // vector from b to a
    double delta[Dim];

    // distance var
    double dist = 0;

    // for all dim
    for (int d = 0; d < Dim; d++)
    {
        // delta in d dim
        delta[d] = a[d] - b[d];

        // squared distance sum
        dist += delta[d] * delta[d];
    }

    // distance square root
    dist =  sqrt(dist);

    // force function
        double f = (exp(- (dist * dist) / force_range ) * force_amount);
       // double f = exp(-dist * dist * 0.01) * 10.0) / dist;


    // compute the force vector for each dim
    for (int d = 0; d < Dim; d++)
    {
        result[d] = delta[d] * f/dist;
    };
}
 */

} // end namespace
