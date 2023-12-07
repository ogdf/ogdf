[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Unreleased

# Porting from Elderberry to current unreleased version

## Requirements
OGDF now requires C++17 features.
We no longer officially support compilers older than gcc 9, clang 9 or Visual Studio 2017 15.8 (MSVC 19.15).

## GraphIO

### SvgPrinter
When their height is equal to their width, `Shape::Triangle` and `Shape::InvTriangle` are now drawn as regular triangles in SVGs (and not as just isosceles ones).

### TikzWriter
`TikzWriter::isCoveredBy()` was removed in favor of `isPointCoveredByNode()` in `geometry.h`.
