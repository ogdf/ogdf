#include <bandit/bandit.h>

#include <ogdf/graphalg/MaxAdjOrdering.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>

using namespace bandit;
using namespace ogdf;

go_bandit([](){
    describe("Maximum Adjacency Orderings ", [](){
        it("should calculate exactly all MAOs.\n", [](){
            for (int N = 4; N <= 9; N++){
                std::cout << "Busy with graphs that have " << N << " nodes." << std::endl;
                Graph P;
                for (int i = 0; i < N; i++){
                    P.newNode();
                }
                MaxAdjOrdering perms;
                ListPure< ListPure< NodeElement *> > allPerms;
                perms.calcAll(&P,&allPerms);

                for (int i = 1; i < 10; i++){
                    Graph G;
                    randomSimpleGraph(G,N,1*N);

                    //make an instance for the MAOs
                    MaxAdjOrdering m;

                    //make structures for saving all MAOs of G
                    ListPure< ListPure< NodeElement *> > MAOs;

                    //calculate them
                    m.calcAll(&G,&MAOs);

                    AssertThat(m.testIfAllMAOs(&G,&MAOs,&allPerms,0),Equals(1));
                }
            }
        });
        it("should calculate MAOs with correct lex-bfs tie breaking.\n", [](){
            for (int N = 10; N <= 20; N++){
                std::cout << "Busy with graphs that have " << N << " nodes." << std::endl;

                for (int i = 1; i < 10; i++){
                    Graph G;
                    randomSimpleGraph(G,N,0.5*N*(N-4));

                    //make an instance for the MAOs
                    MaxAdjOrdering m;

                    //make structures for saving all MAOs of G
                    ListPure< NodeElement *> MAO;

                    //calculate them
                    m.calcBfs(&G,&MAO);

                    AssertThat(m.testIfMAO(&G,&MAO),Equals(1));
                    AssertThat(m.testIfMAOBfs(&G,&MAO),Equals(1));
                }
            }
        });
    });
});
