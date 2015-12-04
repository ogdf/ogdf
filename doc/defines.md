[OGDF](/README.md) » Compiler Defines

# Compiler Defines

Here we list some of the available compiler definitions and macros.

## Defines set by OGDF

| Definition | Location | Description |
|-------------------------------|--------------------|-----------|
| `OGDF_DEBUG`          | build configuration | Perform OGDF assertions.
| `OGDF_DLL`            | build configuration | Building or using ODFG as a DLL.
| `OGDF_INSTALL`        | build configuration | Building ODFG as a DLL.
| `COIN_OSI_CLP`        | build configuration | Clp is the linear program solver (default).
| `COIN_OSI_GRB`        | build configuration | Gurobi is the linear program solver.
| `COIN_OSI_CPX`        | build configuration | CPLEX is the linear program solver.
| `COIN_OSI_CPX`        | build configuration | CPLEX is the linear program solver.
| `OGDF_MEMORY_POOL_TS` | build configuration | OGDF uses the custom thread-safe pool memory manager (default).
| `OGDF_MEMORY_POOL_NTS`    | build configuration | OGDF uses the custom non-thread-safe pool memory manager.
| `OGDF_MEMORY_POOL_MALLOC` | build configuration | OGDF uses the default c++ memory manager.
| `OGDF_SYSTEM_WINDOWS` | `basic/basic.h`  | Compiling for a Windows system.
| `OGDF_SYSTEM_OSX`     | `basic/basic.h`  | Compiling for a Mac OS X system; in this case OGDF_SYSTEM_UNIX is also defined.
| `OGDF_ARCH_X86`       | `basic/System.h` | Compiling for a 32-bit x86 (Intel/AMD) architecture.
| `OGDF_ARCH_X64`       | `basic/System.h` | Compiling for a 64-bit x64 (Intel/AMD) architecture.

# Useful Macros

| Macro | Location | Description |
|-------------------------------|--------------------|-----------|
| `OGDF_ASSERT(expr)`             | `basic/basic.h`      | Assert condition expr in debug builds.
| `OGDF_ASSERT_IF(minLevel,expr)` | `basic/basic.h`      | Assert condition expr in debug builds if debug level is at least minLevel.
| `OGDF_SET_DEBUG_LEVEL(level)`   | `basic/basic.h`      | Set debug level to level.
| `THROW_PARAM(class,param)`      | `basic/exceptions.h` | Throws an exception of type class and passes parameters param; also passes file name and line number of occurrence in source code if OGDF_THROW_WITH_INFO is set.
| `THROW(class)`                  | `basic/exceptions.h` | Throws an exception of type class; also passes file name and line number of occurrence in source code if OGDF_THROW_WITH_INFO is set.
| `OGDF_NEW`                      | `basic/memory.h`     | Should be used to allocated new objects of classes managed by OGDF's memory manager.
| `OGDF_NEW_DELETE`               | `basic/memory.h`     | Adding this macro in a class declaration makes this class managed by OGDF's memory manager.
| `OGDF_MALLOC_NEW_DELETE`        | `basic/memory.h`     | Adding this macro in a class declaration makes this class managed by malloc (as usual), but throws an InsufficientMemoryException when memory allocation failed.
| `OGDF_CHECK_SSE2`               | `basic/System.h`     | Evaluates to true if the current architecture supports SSE2 extensions.
| `OGDF_EXPORT`                   | `basic/basic.h`      | Specifies that a function or class is exported by the OGDF DLL. Set according to the definition of OGDF_INSTALL (OGDF is build as DLL) and OGDF_DLL (OGDF is used as DLL); if none of these is defined (OGDF is build or used as static library), the define expands to an empty string.
| `OGDF_DEPRECATED`               | `basic/basic.h`      | Mark a class / member / function as deprecated.
