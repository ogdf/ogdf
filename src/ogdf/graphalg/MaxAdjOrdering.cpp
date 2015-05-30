/** \file
         * \brief Calculate one or all Maximum Adjacency Ordering(s) of a given simple undirected graph.
         *
         * \author Sebastian Semper
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

#include <ogdf/graphalg/MaxAdjOrdering.h>


namespace ogdf {


MaxAdjOrdering::MaxAdjOrdering(){}
MaxAdjOrdering::~MaxAdjOrdering(){

}

void MaxAdjOrdering::calc(
        const Graph *G,
        ListPure<NodeElement *> *MAO
        ){
    //node count
    int n = G->numberOfNodes();

    //store unsorted nodes
    List<NodeElement *> unsortedNodes;
    G->allNodes(unsortedNodes);

    //neighbourhood counter
    std::vector<int> r(n,0);

    //currently maximal node
    node curMaxNode(*(unsortedNodes.begin()));
    int curMaxVal = 0;

    //last added node
    node lastAdded;

    //add n vertices to M
    for (int i = 0; i < n; i++){

        //add current maximum to end of MAO
        MAO->pushBack(curMaxNode);

        //store the last added
        lastAdded = curMaxNode;

        //delete it from unsorted nodes
        unsortedNodes.del(unsortedNodes.search(lastAdded));

        //set maximal node to first unsorted
        if (i < n-1){
            curMaxNode = unsortedNodes.front();
            curMaxVal = r[curMaxNode->index()];
            for (auto& u : unsortedNodes){
                if (curMaxVal < r[u->index()]){
                    curMaxVal = r[u->index()];
                    curMaxNode = u;
                }
            }
        }

        //edges to iterate over
        edge e;
        forall_adj_edges(e,lastAdded){
            //node at the other side
            node end(e->opposite(lastAdded));

            //search it in unsorted nodes
            ListIterator<NodeElement *> endIt(unsortedNodes.search(end));

            //proceed if unsorted
            if (endIt.valid()){
                //reset neighbourhood
                r[end->index()]++;

                //reset maximal node
                if (r[end->index()] > curMaxVal){
                    curMaxVal = r[end->index()];
                    curMaxNode = end;
                }
            }
        }
    }
}

void MaxAdjOrdering::calcBfs(
        const Graph *G,
        ListPure<NodeElement *> *MAO
        ){
    //node count
    int n = G->numberOfNodes();

    //store unsorted nodes
    ListPure<node> unsortedNodes;
    G->allNodes(unsortedNodes);

    //store tied nodes
    ListPure<node> tiedNodes;
    int curMaxTie = 0;
    node curMaxTieNode(*(unsortedNodes.begin()));

    //neighbourhood counter
    //std::vector<int> r(n,0);
    NodeArray<int> r(*G,0);

    //currently maximal node
    node curMaxNode(*(unsortedNodes.begin()));
    int curMaxVal = 0;

    //last added node
    node lastAdded;

    //add n vertices to M
    for (int i = 0; i < n; i++){

        //add current maximum to end of MAO
        //MAO->pushBack(curMaxNode);
        MAO->pushBack(curMaxTieNode);

        //store the last added
        lastAdded = MAO->back();

        //delete it from unsorted nodes
        unsortedNodes.del(unsortedNodes.search(lastAdded));

        //set maximal node to first unsorted
        if (i < n-1){
            //reset the tied nodes
            tiedNodes.clear();
            tiedNodes.pushBack(unsortedNodes.front());
            curMaxTieNode = unsortedNodes.front();
            curMaxTie = 0;

            //reset the currently maximal node
            curMaxNode = unsortedNodes.front();
            curMaxVal = r[curMaxNode];


            for (auto& u : unsortedNodes){
                //if we find a node tha currently also attains the maximum, we add it to the tiedNodes
                if ((r[u] == curMaxVal)&&(u != tiedNodes.front())){
                    tiedNodes.pushBack(u);
                }
                //if the maximum changes we need to clear the tied nodes
                if (curMaxVal < r[u]){
                    tiedNodes.clear();
                    tiedNodes.pushBack(u);
                    curMaxVal = r[u];
                    curMaxNode = u;
                }
            }
            //now calc the lex-bfs-value for every tied node
            for (auto& t : tiedNodes){
                edge e;
                int tieVal = 0;
                forall_adj_edges(e,t){
                    ListIterator<node> opIt = MAO->search(e->opposite(t));
                    if (opIt.valid()){
                        tieVal += pow(2,MAO->size()-MAO->pos(opIt));
                    }
                }
                //update the currently maximum tied node
                if (tieVal > curMaxTie){
                    curMaxTieNode = t;
                    curMaxTie = tieVal;
                }
            }
        }

        //edges to iterate over
        edge e;
        forall_adj_edges(e,lastAdded){
            //node at the other side
            node end(e->opposite(lastAdded));

            //search it in unsorted nodes
            ListIterator<NodeElement *> endIt(unsortedNodes.search(end));

            //proceed if unsorted
            if (endIt.valid()){
                //reset neighbourhood
                r[end]++;
                bool tieUpdated = 0;
                if (r[end] == curMaxVal){
                    tiedNodes.pushBack(end);
                    tieUpdated = 1;
                }
                //reset maximal node
                if (r[end] > curMaxVal){
                    tiedNodes.clear();
                    tiedNodes.pushBack(end);
                    curMaxVal = r[end];
                    curMaxNode = end;
                    tieUpdated = 1;
                }

                //if anything changed we need to recalc all tied nodes
                //TODO only do what is neccesary
                if (tieUpdated){
                    curMaxTieNode = end;
                    curMaxTie = 0;
                    for (auto& t : tiedNodes){
                        edge e2;
                        int tieVal = 0;
                        forall_adj_edges(e2,t){
                            ListIterator<node> opIt = MAO->search(e2->opposite(t));
                            if (opIt.valid()){
                                tieVal += pow(2,MAO->size()-MAO->pos(opIt));
                            }
                        }
                        //update the currently maximum tied node
                        if (tieVal > curMaxTie){
                            curMaxTieNode = t;
                            curMaxTie = tieVal;
                        }
                    }
                }
            }
        }

        std::cout << "Tied nodes with maximal tie value -" << curMaxTie << "- ang the nodes: ";
        for (auto& t : tiedNodes){
            std::cout << t->index() << ",";
        }
        std::cout <<" and  node "<< curMaxTieNode->index() << " wins." << std::endl;

    }
}

void MaxAdjOrdering::calc(
        const Graph *G,
        ListPure<NodeElement *> *MAO,
        ListPure< ListPure<EdgeElement *> > *Forests
        ){
    //node count
    int n(G->numberOfNodes());

    //store unsorted nodes
    List<NodeElement *> unsortedNodes;
    G->allNodes(unsortedNodes);

    //neighbourhood counter
    std::vector<int> r(n,0);

    //currently maximal node
    node curMaxNode(*(unsortedNodes.begin()));
    int curMaxVal(0);

    //last added node
    node lastAdded;

    //add n vertices to M
    for (int i = 0; i < n; i++){

        //add current maximum to end of MAO
        MAO->pushBack(curMaxNode);

        //store the last added
        lastAdded = curMaxNode;

        //delete it from unsorted nodes
        unsortedNodes.del(unsortedNodes.search(lastAdded));

        //set maximal node to currently maximal unsorted
        if (i < n-1){
            curMaxNode = unsortedNodes.front();
            curMaxVal = r[curMaxNode->index()];
            for (auto& u : unsortedNodes){
                if (curMaxVal < r[u->index()]){
                    curMaxVal = r[u->index()];
                    curMaxNode = u;
                }
            }
        }

        //edges to iterate over
        edge e;
        forall_adj_edges(e,lastAdded){
            //node at the other side
            node end(e->opposite(lastAdded));

            //search it in unsorted nodes
            ListIterator<NodeElement *> endIt(unsortedNodes.search(end));

            //proceed if unsorted
            if (endIt.valid()){
                //correct neighbourhood
                int r_(++r[end->index()]);

                if (r_ > curMaxVal){
                    curMaxVal = r_;
                    curMaxNode = end;
                }
                if (r_ >= Forests->size()){
                    Forests->pushBack(ListPure<EdgeElement *>());
                }
                (*(Forests->get(r_-1))).pushBack(e);
            }
        }
    }
}

void MaxAdjOrdering::calc(
        const Graph *G,
        NodeElement *s,
        ListPure<NodeElement *> *MAO
        ){
    //node count
    int n = G->numberOfNodes();

    //store unsorted nodes
    List<NodeElement *> unsortedNodes;
    G->allNodes(unsortedNodes);

    //neighbourhood counter
    std::vector<int> r(n,0);

    //currently maximal node
    node curMaxNode = s;
    int curMaxVal = 0;

    //last added node
    node lastAdded;

    //add n vertices to M
    for (int i = 0; i < n; i++){

        //add current maximum to end of MAO
        MAO->pushBack(curMaxNode);

        //store the last added
        lastAdded = curMaxNode;

        //delete it from unsorted nodes
        unsortedNodes.del(unsortedNodes.search(lastAdded));

        if (i < n-1){
            curMaxNode = unsortedNodes.front();
            curMaxVal = r[curMaxNode->index()];
            for (auto& u : unsortedNodes){
                if (curMaxVal < r[u->index()]){
                    curMaxVal = r[u->index()];
                    curMaxNode = u;
                }
            }
        }

        //edges to iterate over
        edge e;
        forall_adj_edges(e,lastAdded){
            //node at the other side
            node end(e->opposite(lastAdded));

            //search it in unsorted nodes
            ListIterator<NodeElement *> endIt(unsortedNodes.search(end));

            //proceed if unsorted
            if (endIt.valid()){
                //correct neighbourhood
                r[end->index()]++;

                if (r[end->index()] > curMaxVal){
                    curMaxVal = r[end->index()];
                    curMaxNode = end;
                }
            }
        }
    }
}

void MaxAdjOrdering::calc(
        const Graph *G,
        NodeElement *s,
        ListPure<NodeElement *> *MAO,
        ListPure< ListPure<EdgeElement *> > *Forests
        ){
    //node count
    int n = G->numberOfNodes();

    //store unsorted nodes
    List<NodeElement *> unsortedNodes;
    G->allNodes(unsortedNodes);

    //neighbourhood counter
    std::vector<int> r(n,0);

    //currently maximal node
    node curMaxNode = s;
    int curMaxVal = 0;

    //last added node
    node lastAdded;

    //add n vertices to M
    for (int i = 0; i < n; i++){

        //add current maximum to end of MAO
        MAO->pushBack(curMaxNode);

        //store the last added
        lastAdded = curMaxNode;

        //delete it from unsorted nodes
        unsortedNodes.del(unsortedNodes.search(lastAdded));

        //set maximal node to first unsorted
        if (i < n-1){
            curMaxNode = unsortedNodes.front();
            curMaxVal = r[curMaxNode->index()];
            for (auto& u : unsortedNodes){
                if (curMaxVal < r[u->index()]){
                    curMaxVal = r[u->index()];
                    curMaxNode = u;
                }
            }
        }

        //edges to iterate over
        edge e;
        forall_adj_edges(e,lastAdded){
            //node at the other side
            node end(e->opposite(lastAdded));

            //search it in unsorted nodes
            ListIterator<NodeElement *> endIt(unsortedNodes.search(end));

            //proceed if unsorted
            if (endIt.valid()){
                //correct neighbourhood
                int r_(++r[end->index()]);

                if (r[end->index()] > curMaxVal){
                    curMaxVal = r[end->index()];
                    curMaxNode = end;
                }

                if (r_ >= Forests->size()){
                    Forests->pushBack(ListPure<EdgeElement *>());
                }
                (*(Forests->get(r_-1))).pushBack(e);
            }
        }
    }
}

void MaxAdjOrdering::calcAll(
        const Graph *G,
        ListPure< ListPure<NodeElement *> > *MAOs
        ){
    //initialize backtrackstack
    ListPure<NodeElement *> nodes;
    G->allNodes(nodes);

    //first step in recursion. every node is an option for the first one
    //in the ordering. so start the recursion once for every node
    for(NodeElement *it : nodes){
        ListPure<NodeElement *> start;
        ListPure<NodeElement *> unsorted = nodes;
        unsorted.del(unsorted.search(it));
        start.pushBack(it);
        m_calcAllMAOs_recursion(
                    G->numberOfNodes(),
                    start,
                    unsorted,
                    std::vector<int>(G->numberOfNodes(),0),
                    MAOs
                    );
    }
}


void MaxAdjOrdering::m_calcAllMAOs_recursion(
        int n,
        ListPure<NodeElement *> currentOrder,
        ListPure<NodeElement *> currentUnsorted,
        std::vector<int> r,
        ListPure< ListPure<NodeElement *> > *MAOs
        ){
    if (currentUnsorted.empty()){
        //one MAO is done!
        MAOs->pushBack(currentOrder);

        //go back up into recursion
        return;
    }

    //store the last node in current order
    node lastAdded = currentOrder.back();

    //if we want all maos, we have to store ALL possible next nodes
    ListPure<NodeElement *> maxNodes;

    //choose the first maxValue as the first value in the unsorted
    int maxValue(
                r[(currentUnsorted.front())->index()]
            );

    for (auto& u : currentUnsorted){
        if (maxValue < r[u->index()]){
            maxValue = r[u->index()];
        }
    }

    //add all nodes that have this value
    for(NodeElement * it : currentUnsorted){
        if (r[it->index()] == maxValue){
            maxNodes.pushBack(it);
        }
    }

    //edges to iterate over
    edge e;
    forall_adj_edges(e,lastAdded){
        //node at the other side
        node end(e->opposite(lastAdded));
        ListIterator<NodeElement *> endIt(currentUnsorted.search(end));

        //if is unsorted
        if (endIt.valid()){
            int endInd((*endIt)->index());

            //increase value of neighborhood
            r[endInd]++;

            //if it is the current maximum, add it to the list
            if (r[endInd] == maxValue){
                maxNodes.pushBack(end);
            }

            //if it is larger than the current maximum
            if (r[endInd] > maxValue){

                //reset maximum value
                maxValue = r[endInd];

                //clear maximum list
                maxNodes.clear();

                //add the current node
                maxNodes.pushBack(end);
            }
        }
    }

    //go deeper into recursion for every possible node in maxNodes
    for(NodeElement *it : maxNodes){
        ListPure<NodeElement *> nextOrder = currentOrder;
        ListPure<NodeElement *> nextUnsorted = currentUnsorted;

        //the current node is the next one in the next calculated order
        nextOrder.pushBack(it);

        //the current node needs to be removed from the unsorted for the next step
        nextUnsorted.del(nextUnsorted.search(it));
        m_calcAllMAOs_recursion(
                    n,
                    nextOrder,
                    nextUnsorted,
                    r,
                    MAOs
                    );
    }
}

void MaxAdjOrdering::calcAll(
        const Graph *G,
        ListPure< ListPure<NodeElement *> > *MAOs,
        ListPure< ListPure< ListPure <EdgeElement *> > > *Fs
        ){
    //initialize backtrackstack
    ListPure<NodeElement *> nodes;
    G->allNodes(nodes);

    //first step in recursion. every node is an option for the first one
    //in the ordering. so start the recursion once for every node
    //but just leave the forests empty at first
    for(NodeElement *it : nodes){
        ListPure<NodeElement *> start;
        ListPure<NodeElement *> unsorted = nodes;
        unsorted.del(unsorted.search(it));
        start.pushBack(it);
        m_calcAllMAOs_recursion(
                    G->numberOfNodes(),
                    start,
                    ListPure< ListPure <EdgeElement *> >(),
                    unsorted,
                    std::vector<int>(G->numberOfNodes(),0),
                    MAOs,
                    Fs
                    );
    }
}

void MaxAdjOrdering::m_calcAllMAOs_recursion(
        int n,
        ListPure<NodeElement *> currentOrder,
        ListPure< ListPure <EdgeElement *> > currentForest,
        ListPure<NodeElement *> currentUnsorted,
        std::vector<int> r,
        ListPure< ListPure<NodeElement *> > *MAOs,
        ListPure< ListPure< ListPure <EdgeElement *> > > *Fs
        ){
    if (currentUnsorted.empty()){
        //one MAO is done!
        MAOs->pushBack(currentOrder);
        Fs->pushBack(currentForest);
        //go back up into recursion
        return;
    }

    //store the last node in current order
    node lastAdded = currentOrder.back();

    //if we want all maos, we have to store ALL possible next nodes
    ListPure<NodeElement *> maxNodes;

    //choose the first maxValue as the first value in the unsorted
    int maxValue(
                r[(currentUnsorted.front())->index()]
            );
    for (auto& u : currentUnsorted){
        if (maxValue < r[u->index()]){
            maxValue = r[u->index()];
        }
    }

    //add all nodes that have this value
    for(NodeElement * it : currentUnsorted){
        if (r[it->index()] == maxValue){
            maxNodes.pushBack(it);
        }
    }

    //edges to iterate over
    edge e;
    forall_adj_edges(e,lastAdded){
        //node at the other side
        node end(e->opposite(lastAdded));
        ListIterator<NodeElement *> endIt(currentUnsorted.search(end));

        //if is unsorted
        if (endIt.valid()){
            //increase value of neighborhood and store it
            float r_(++r[((*endIt)->index())]);

            //if it is the current maximum, add it to the list
            if (r_ == maxValue){
                maxNodes.pushBack(end);
            }

            //if it is larger than the current maximum
            if (r_ > maxValue){

                //reset maximum value
                maxValue = r_;

                //clear maximum list
                maxNodes.clear();

                //add the current node
                maxNodes.pushBack(end);
            }

            //depending on the last node - populate the forest accordingly
            if (r_ >= currentForest.size()){
                currentForest.pushBack( ListPure <EdgeElement *> ());
            }
            (*currentForest.get(r_-1)).pushBack(e);
        }
    }

    //go deeper into recursion for every possible node in maxNodes
    for(NodeElement *it : maxNodes){
        ListPure<NodeElement *> nextOrder = currentOrder;
        ListPure<NodeElement *> nextUnsorted = currentUnsorted;

        //the current node needs to be removed from the unsorted for the next step
        nextUnsorted.del(nextUnsorted.search(it));

        //the current node is the next one in the next calculated order
        nextOrder.pushBack(it);

        m_calcAllMAOs_recursion(
                    n,
                    nextOrder,
                    currentForest,
                    nextUnsorted,
                    r,
                    MAOs,
                    Fs
                    );
    }
}

bool MaxAdjOrdering::testIfMAO(const Graph *G, ListPure<NodeElement *> *Ordering)
{
    unsigned int i = 0;
    unsigned int n = Ordering->size();
    std::vector<unsigned int> r(Ordering->size(),0);
    edge e;
    node op;
    ListPure<NodeElement *> tested;
    for (auto& o : *Ordering){
        tested.pushBack(o);
        forall_adj_edges(e,o){
            op = e->opposite(o);
            //check if edge goes to the right
            if (!tested.search(op).valid()){
                //increment the neighbourhood counter
                r[op->index()]++;
            }
        }
        if (i < n-1){
            /**go through all following nodes and check if
                                                         * neighbourhood is bigger than then one in the
                                                         * ordering. If yes - return false because no MAO.
                                                         */
            for (ListIterator<NodeElement* > next = Ordering->get(i+1); next != Ordering->end(); next++){
                if (r[(*next)->index()] > r[(*(Ordering->get(i+1)))->index()]){
                    return 0;
                }
            }
        }
        i++;
    }
    return 1;
}

bool MaxAdjOrdering::testIfAllMAOs(
        const Graph *G,
        ListPure< ListPure< NodeElement *> > *Orderings,
        ListPure< ListPure< NodeElement *> > *Perms,
        bool verbose
        ){
    ListPure<NodeElement *> nodes;
    G->allNodes(nodes);
    int n = nodes.size();
    ListPure<NodeElement *> testPerm;

    for (auto& p : *Perms){
        testPerm.clear();
        //generate nodelist of G from permutation
        for (int i = 0; i < n; i++){
            int index = (*(p.get(i)))->index();
            if(verbose){
                std::cout << index+1;
            }
            testPerm.pushBack(
                        *(nodes.get(index))
                        );
        }

        //check if testPerm is a MAO - this way we find all MAOs
        if (testIfMAO(G,&testPerm)){
            //if we don't find the testPerm in our provided list, we will not have generated
            //all MAOs
            ListIterator<ListPure <NodeElement *>> pIt = Orderings->search(testPerm);
            if (!pIt.valid()){
                if(verbose){
                    std::cout << " - [ERR] Found permutation that is MAO and not list!" << std::endl;
                }
                return 0;
            }
            else {
                if(verbose){
                    std::cout << " - [INF] Found permutation that is MAO and in list!" << std::endl;
                }
            }
        }
        else {
            //if we find the testPerm in the list we did calculate to many MAOs
            ListIterator<ListPure <NodeElement *>> pIt = Orderings->search(testPerm);
            if (!pIt.valid()){
                if(verbose){
                    std::cout << " - [INF] Found permutation that no MAO and not list!" << std::endl;
                }
            }
            else {
                if(verbose){
                    std::cout << " - [ERR] Found permutation that is no MAO and in list!" << std::endl;
                }
                return 0;
            }
        }
    }
    if(verbose){
        std::cout << std::endl;
    }

    return 1;
}

void MaxAdjOrdering::visualize(
        GraphAttributes *GA,
        ListPure<NodeElement *> *MAO
        ){
    const Graph& G = GA->constGraph();
    ogdf::List<ogdf::NodeElement *> nodes;
    G.allNodes(nodes);

    ogdf::LinearLayout l(600,*MAO);
    l.setCustomOrder(1);
    l.call(*GA);

    int k = 1;
    for(auto& n : *MAO){
        GA->height(n) = 15;
        GA->width(n) = 15;
        GA->label(n) = std::to_string(k++);
        GA->shape(n) = shEllipse;
        GA->strokeColor(n) = Color(Color::Black);
        GA->fillColor(n) = Color(Color::Red);
    }
}

void MaxAdjOrdering::visualize(
        GraphAttributes *GA,
        ListPure<NodeElement *> *MAO,
        ogdf::ListPure< ogdf::ListPure <ogdf::EdgeElement *> > *F
        ){
    const Graph& G = GA->constGraph();
    ogdf::List<ogdf::NodeElement *> nodes;
    G.allNodes(nodes);

    ogdf::LinearLayout l(140*nodes.size(),*MAO);
    l.setCustomOrder(1);
    l.call(*GA);

    int k = 1;
    for(auto& n : *MAO){
        GA->height(n) = 15;
        GA->width(n) = 30;
        GA->label(n) = std::to_string(k++) + std::string(",") + std::to_string(n->index()+1);
        GA->shape(n) = shEllipse;
        GA->strokeColor(n) = Color(Color::Black);
        GA->fillColor(n) = Color(Color::Red);
    }

    k = 1;
    for (auto& f : *F){
        for (auto& e : f){
            GA->strokeWidth(e) = 2*k;
            GA->arrowType(e) = eaNone;
        }
        k++;
    }
}
} // end namespace ogdf

