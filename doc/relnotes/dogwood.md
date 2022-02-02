[OGDF](../../README.md) » [Release Notes](../relnotes.md) » Dogwood

# OGDF Dogwood (2022.02)

Released 2022-02-02.

As always, this release includes many new features, bug fixes and documentation improvements.
The most important changes are listed below.
Note that we are dropping the official support of some older compilers and are aiming to make use of C++17 features in the future.
For more details, take a look at the [porting guide](../porting/dogwood.md).

Noteworthy changes:
* compilation:
    * fixes for compilation with gcc 11, clang 11, AppleClang 12, and ClangCL on MSVC
    * fixes for compilation on ARM architectures
    * installation of shared libraries now possible
    * new cmake option `OGDF_FULL_DOC` to include {`test/`,`src/`} doc in doxygen output
    * new cmake option `OGDF_ENABLE_CLANG_TIDY` to enable static analysis via clang-tidy
    * `OGDF_` prepended to all commented out OGDF-macros (in case you enable them...)
* basic graph functionality:
    * new `GraphObjectContainer<E>::permute()`, allows e.g. `G.nodes.permute()`
    * new `GraphList::empty()`
    * new optional parameter for `contract()` to keep self-loops
    * `splitNode(a, a)` is handled correctly
    * `moveAdj()` handles self-loops correctly
* embeddings and dual graph:
    * new `operator<<()` for `face`s
    * `CombinatorialEmbedding::splitFace()` can create self-loops
    * `CombinatorialEmbedding::joinFaces()` accepts a bridge as an argument
    * `EmbedderMaxFace` now works correctly
    * new `DynamicDualGraph` allowing dynamic changes of a dual graph
* `GraphCopy`:
    * `GraphCopy::removeEdgePathEmbedded()` now removes degree-1-nodes
    * new overloads of `GraphCopy::{insert,remove}EdgePathEmbedded()` can change a `DynamicDualGraph`
    * new `GraphCopy` methods to detect and remove non-simple crossings
    * `GraphCopy[Simple]::init()` now actually clears the graph first
    * empty `GraphCopySimple` can be constructed, similar to `GraphCopy`
    * `GraphCopy` of an empty `Graph` is now initialized correctly
* `GraphAttributes`:
    * `idNode()` now returns `v->index()` if no user id is set
* `CoinManager`:
    * `updateLogging()` replaces `logging()`, propagates log level to COIN
* (graph) algorithms:
    * new `findCutVertices()` to extract the cut vertices of a graph
    * `triangulate()` now runs in linear time (previously: quadratic)
    * `Dijkstra`: new optional parameters for early termination
    * `STNumbering` no longer uses recursion, large instances do not cause a stack overflow
    * new `BCTree::initNotConnected()` for a given subset of graph nodes/components
    * new `DisjointSets::init()` allows reinitialization
    * new crossing minimization heuristics employing star insertion:
        * `PlanarizerStarReinsertion`
        * `PlanarizerMixedInsertion`
        * `PlanarizerChordlessCycle`
    * `SubgraphPlanarizer` removes non-simple crossings
* `GraphIO`:
    * `GraphIO::read()` now determines format via file extension
        * `GraphIO::read(G, "input.graphml", GraphIO::read)` for old behavior
    * new `GraphIO::readTsplibXml()` for reading tsplib instances in XML format
    * `drawSVG()` now draws cluster labels
    * `drawSVG()` now draws arrows heads correctly
    * `readDOT()` now recognizes numeric literals correctly
    * `readGML()`/`writeGML()` now respects `idNode()` attribute
    * various fixes for `ClusterGraph` writing
* layouts:
    * `FFPLayout` and `ComponentSplitterLayout` respect a given embedding if possible
    * `FMMMLayout`: high-level options no longer reset low-level options
        * new `FMMMLayout::resetOptions()` allows for manual option reset
    * `HierarchyLayoutModule` no longer resets node attributes
    * `PivotMDS` computes 3D coordinates when `GraphAttributes::threeD` is set
    * `SugiyamaLayout` no longer reverses bend points
    * `SugiyamaLayout` can handle `ClusteredGraphAttributes`

This release contains (big and small) contributions by Antoine Lambert, Finn Stutzenstein,
Hendrik Brückler, Ivo Hedtke, Jöran Schierbaum, Mario Emmenlauer, Matthias Pfretzschner,
Max Ilsen, Simon Dominik "Niko" Fink, Stephan Beyer, Thomas Klein, Thomas Roehr,
Vadim Zabermakh, neotechllc and xuanjueheshang on GitHub. Thanks a lot to all contributors!
