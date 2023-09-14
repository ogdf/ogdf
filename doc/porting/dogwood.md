[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Dogwood

# Porting from Catalpa to Dogwood

## Requirements

We no longer officially support compilers older than gcc 7, clang 6 or Visual Studio 2017 15.8 (MSVC 19.15). Further, we require CMake version 3.8 or newer.
This allows us to ignore certain quirks of these older compilers and even make use of C++17 features in the future.
Some parts of the OGDF may still be usable even when using older compiler versions, but we will not pay any extra effort to keep it that way.

## API Changes

### ArrayBuffer
`ArrayBuffer::compactMemcpy()` is removed in favor of `ArrayBuffer::compactCopy()`, which now automatically uses `memcpy` if the array elements are trivially copiable.

### CoinManager logging
`CoinManager::logging(OsiSolverInterface* osi, bool logMe)` is removed in favor of `updateLogging(OsiSolverInterface* osi)`,
which propagates changes made to the newly introduced `ogdf::Logger CoinManager::CoinLog` to the COIN logging facilities
in a more fine-grained way.

### Default value for GraphAttributes::idNode
When using `GraphAttributes` with `nodeId` enabled, the `GraphAttributes::idNode()` method returns a given node's current
index (rather than the old default value of `-1`) if no user-specified index has been explicitly set.

### File-Type guessing
Using `read(G, "input.graphml")` on a file that actually contains data in another format no longer works,
as we no longer try to auto-detect the file type if we can deduce it from the file's extension.
If you want the old behavior ignoring the file's extension, use `read(G, "input.graphml", GraphIO::read)`
to always try all known formats until one is found that successfully parses the file.

### FMMM options
`fmmm.useHighLevelOptions(true)` no longer resets all low-level options (see the method's documentation for more information).
If you want to emulate the previous behavior, simply use `fmmm.resetOptions(); fmmm.useHighLevelOptions(true);`.
