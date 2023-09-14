[OGDF](../../README.md) » [Release Notes](../relnotes.md) » Elderberry

# OGDF Elderberry (2023.09)

Released 2023-09-14.

This release marks the beginning of public OGDF development on GitHub! Hooray!
Finally, bug fixes and new features will be directly available to the public.
While certain research branches will still have to be kept in a private
repository until the corresponding scientific results have been published, the
bulk of the development now happens on GitHub.
Most importantly, this means that members of the GitHub community can
participate in the development of the OGDF in a more direct manner.

We also want to highlight [ogdf-python](https://github.com/ogdf/ogdf-python):
These python bindings for the OGDF have proved to work quite well, also in
conjunction with Jupyter Notebook. Try them out!
Further, we now have [conan](https://conan.io/center/recipes/ogdf) and
[vcpkg](https://vcpkg.roundtrip.dev/ports/ogdf) packages for the OGDF available.

## Major Additions Regarding Public OGDF Development

The following new additions will facilitate OGDF development:
* new [dev-guide](../dev-guide.md), inviting you to contribute!
* code style via clang-format (see `.clang-format` in the root directory):
    * new script `test_clang_format.sh` for easy code style enforcement, using
      a public docker image when the correct version is not installed locally
    * accordingly reformatted code base
* continuous integration (CI) using GitHub Actions:
    * [publicly available docker images](https://hub.docker.com/u/ogdf) for the
       CI, and scripts to build these images
    * scripts for compiling and running OGDF
        * for Linux, Mac and Windows (the former two using ccache)
        * for different build types (Release, Debug)
        * for different compilers (gcc, clang)
        * with or without Gurobi as the ILP solver
    * scripts for style checks beyond clang-format

## New Features and Bug Fixes

Concerning the main code base, this release includes many new features, bug
fixes, tests and documentation improvements. Most noteworthy are:
* compilation fixes, including:
    * C++20 compatibility
    * ARM64 Windows builds
    * cross-compilation for Windows on Linux
* new geometric crossing minimization functionality:
    * `CrossingMinimalPositionFast`
    * `GeometricEdgeInsertion`, `GeometricVertexInsertion`, and `VertexMovement`
    * requires compilation with [CGAL](https://www.cgal.org/), for more info
      consult the corresponding [README](../../include/ogdf/geometric/README.md)
* new `MaximumDensitySubgraph`
* new minimum cut algorithms:
    * `MinimumCutModule`
    * `MinimumCutNagamochiIbaraki`
    * `MinimumCutStoerWagner`, a rewritten `MinCut` with performance in mind
* new planar separator algorithms:
    * `PlanarSeparatorModule`
    * `SeparatorDual`
    * `SeparatorDualFC`
    * `SeparatorHarPeled`
    * `SeparatorLiptonTarjan`
    * `SeparatorLiptonTarjanFC`
* new spanner algorithms:
    * `SpannerModule`
    * `SpannerBasicGreedy`
    * `SpannerBaswanaSen`
    * `SpannerBerman` and `SpannerBermanDisconnected`
    * `SpannerElkinNeiman`
    * `SpannerIteratedWrapper`
    * `SpannerKortsarzPeleg`
* new `findMaximumCardinalityMatching()`
* new `MaxAdjOrdering::calcForest()` and `MaxAdjOrdering::visualize()`
* new `globeGraph()` generator
* `GraphIO`:
    * new `drawTikz()` allowing to write out graphs as LaTeX/TikZ
    * new `GraphAttributes::isUniform()` to test attribute uniformity
    * `drawSVG()` no longer dashes arrow heads of dashed arrows
    * `drawSVG()` now correctly uses its given settings
    * `readDL()` no longer fails when reading nodes with index n
    * `readGML()` no longer crashes on malformed index data
    * `readGraphML()` now accepts yFiles yEd files
* layouts:
    * rewritten `RadialTreeLayout` now ensures planarity
    * new `setForcing2DLayout()` for `PivotMDS` and `StressMinimization`,
      ensuring a 2D-layout even when `GraphAttributes::threeD` is set

## Contributors

This release contains (big and small) contributions by Dominik Potulski,
Finn Stutzenstein, Hendrik Brückler, Ivo Hedtke, Jens Schmidt, Jöran Schierbaum,
Kassian Köck, Marcel Radermacher, Mark Scheibner, Max Ilsen, Mirko H. Wagner,
Sascha Alder, Simon Dominik "Niko" Fink, Stephan Beyer, Thomas Klein, and
modmuss50 on GitHub.
