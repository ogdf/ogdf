/** \file
 * \brief Tests for the basic disjoint sets class
 *
 * \author Mirko Wagner
 *
 * \par License:
 * This file is part of the Open myGraPath Halving, Drawing Framework (OGDF).
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
 * Software Foundation, Inc., 51 FRank,lin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <bandit/bandit.h>
#include <ogdf/basic/DisjointSets.h>
#include <resources.h>
#include <ogdf/basic/exceptions.h>

using namespace ogdf;
using namespace bandit;

template<typename DisjointSetsClass>
static void registerTestSuite(const string typeName)
{
    describe(typeName,[&](){
        DisjointSetsClass *disjointSets;
        int sets[42];

        before_each([&](){
            disjointSets = new DisjointSetsClass();
             for(int i = 0; i < 42; i++){
                sets[i] = disjointSets->makeSet();
            }
        });

        after_each([&](){
            delete disjointSets;
        });

        it("assigns valid set id's", [&](){
            for(int i : sets){
                AssertThat(i, IsGreaterThan(-1));
            }
        });

        it("is initialized", [&](){
            DisjointSetsClass *emptydisjointSets;
            emptydisjointSets= new DisjointSetsClass();
            AssertThat(emptydisjointSets->getNumberOfElements(), Equals(0));
            AssertThat(emptydisjointSets->getNumberOfSets(), Equals(0));

        });

        it("can be filled", [&](){
            AssertThat(disjointSets->getNumberOfElements(), Equals(42));
            AssertThat(disjointSets->getNumberOfSets(), Equals(42));
        });

        it("unifies two disjoint sets and doesn't unify two joined sets", [&](){
            AssertThat(disjointSets->quickUnion(sets[2], sets[1]), IsTrue());
            AssertThat(disjointSets->quickUnion(sets[0], sets[2]), IsTrue());
            AssertThat(disjointSets->quickUnion(sets[0], sets[1]), IsFalse());
        });

        it("returns the same id for every item of a unified superset", [&](){
            AssertThat(disjointSets->getRepresentative(sets[13]), Equals(sets[13]));
            AssertThat(disjointSets->getRepresentative(sets[13]), Equals(disjointSets->find(sets[13])));
            disjointSets->quickUnion(sets[1], sets[2]);
            disjointSets->quickUnion(sets[2], sets[3]);
            disjointSets->quickUnion(sets[1], sets[4]);
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[2])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[3])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[4])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[1])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[2])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[3])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[4])));
            AssertThat(sets[5], !Equals(disjointSets->find(sets[4])));
            AssertThat(sets[5], !Equals(disjointSets->getRepresentative(sets[4])));
        });

        it("returns the same id for every item of a linked superset", [&](){
            AssertThat(disjointSets->getRepresentative(13), Equals(13) && Equals(disjointSets->find(13)));
            disjointSets->link(sets[1], sets[2]);
            disjointSets->link(disjointSets->find(sets[2]), sets[3]);
            disjointSets->link(disjointSets->find(sets[1]), sets[4]);
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[2])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[3])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->getRepresentative(sets[4])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[1])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[2])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[3])));
            AssertThat(disjointSets->getRepresentative(sets[1]), Equals(disjointSets->find(sets[4])));
            AssertThat(sets[5], !Equals(disjointSets->find(sets[4])));
            AssertThat(sets[5], !Equals(disjointSets->getRepresentative(sets[4])));
        });

        it("tracks the number of elements", [&](){
            AssertThat(disjointSets->getNumberOfElements(), Equals(42));
            disjointSets->quickUnion(sets[1], sets[2]);
            disjointSets->quickUnion(sets[1], sets[2]);
            disjointSets->link(disjointSets->getRepresentative(sets[1]), disjointSets->find(sets[3]));
            disjointSets->link(disjointSets->find(sets[2]), disjointSets->getRepresentative(sets[3]));
            AssertThat(disjointSets->getNumberOfElements(), Equals(42));
            disjointSets->makeSet();
            AssertThat(disjointSets->getNumberOfElements(), Equals(43));
        });

        it("tracks the number of sets when using link", [&](){
            int successfullLinkCounter = 0;
            successfullLinkCounter += disjointSets->link(sets[1], sets[2]) != -1;
            successfullLinkCounter += disjointSets->link(disjointSets->find(sets[2]), sets[3]) != -1;
            successfullLinkCounter += disjointSets->link(disjointSets->find(sets[1]), sets[4]) != -1;
            successfullLinkCounter += disjointSets->link(disjointSets->find(sets[4]), disjointSets->find(sets[3])) != -1;
            for(int i = 0; i < 42; i++){
                successfullLinkCounter += disjointSets->link(disjointSets->find(sets[4]), disjointSets->find(sets[3])) != -1;
            }
            AssertThat(disjointSets->getNumberOfSets(), Equals(42 - successfullLinkCounter));
            AssertThat(successfullLinkCounter, IsLessThan(42));
        });

        it("tracks the number of sets when using quickUnion", [&](){
            int successfullUnionCounter = 0;
            successfullUnionCounter += disjointSets->quickUnion(sets[1], sets[2]);
            successfullUnionCounter += disjointSets->quickUnion(sets[2], sets[3]);
            successfullUnionCounter += disjointSets->quickUnion(sets[1], sets[4]);
            for(int i=0;i<42;i++){
                successfullUnionCounter += disjointSets->quickUnion(sets[4], sets[3]);
            }

            AssertThat(disjointSets->getNumberOfSets(), Equals(42-successfullUnionCounter));
        });

#ifdef OGDF_DEBUG
        it("throws an exception, if the user tries to link two non-maximal disjoint sets", [&](){
            disjointSets->link(sets[3], sets[4]);
            int notMaximalSet = ( disjointSets->getRepresentative(sets[3]) == sets[4] ? sets[3] : sets[4] );
            AssertThrows(PreconditionViolatedException, disjointSets->link(notMaximalSet, sets[5]));
        });

        it("detects invalid set ids", [&](){
            AssertThrows(PreconditionViolatedException, disjointSets->find(-1));
            AssertThrows(PreconditionViolatedException, disjointSets->getRepresentative(-1));
            int notASetId=0;
            int max=0;
            for(int i : sets){
                max = ( i > max ? i : max );
            }
            notASetId = max+1;
            AssertThrows(PreconditionViolatedException, disjointSets->find(notASetId));
            AssertThrows(PreconditionViolatedException, disjointSets->getRepresentative(notASetId));
        });
#endif
	});
}

go_bandit([](){
describe("Disjoint Sets", [](){
    registerTestSuite<DisjointSets<>>("Default");
    registerTestSuite<DisjointSets<LI, PC, Rem>>("Linking by Index, Path Compression, Rem's Algorithm");
    registerTestSuite<DisjointSets<LR, PS, TvL>>("Linking by Rank, Path Splitting, Tarjan and van Leeuwen's Algorithm");
    registerTestSuite<DisjointSets<NL, R1, IR0>>("No Linking, Reversal Type 1, Interleaved Reversal Type 0");
    registerTestSuite<DisjointSets<LI, PH, IPSPC>>("Linking by Index, Path Halving, Interleaved Path Splitting Path Compression");
    registerTestSuite<DisjointSets<LS, CO, NI>>("Linking by Size, Collapsing, No Interleaving");
    registerTestSuite<DisjointSets<NL, NF, NI>>("No Linking, No Compression, No Interleaving");
});
});
