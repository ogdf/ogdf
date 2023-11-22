[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Unreleased

# Porting from Elderberry to current unreleased version

## GraphIO

### SvgPrinter
When their height is equal to their width, `Shape::Triangle` and `Shape::InvTriangle` are now drawn as regular triangles in SVGs (and not as just isosceles ones).

### TikzWriter
`TikzWriter::isCoveredBy()` was removed in favor of `isPointCoveredByNode()` in `geometry.h`.
