[OGDF](../README.md) » [Developer's Guide](dev-guide.md) » Coding Standards

# Coding Standards {#standards}

This document is about coding standards for the OGDF.
Every OGDF developer is required to read and internalize these coding standards.

Note that the OGDF grew for a long time without or with different and looser coding
standards, so there are a lot of places in the code where these standards
are not met.

----

## Best Practices

The main objective is to write good code.
A terse summary of what good code is, can be found in
[a blog post by Jesse Hill](https://spin.atomicobject.com/2015/03/04/good-code/).
In [The Pragmatic Programmer Quick Reference Guide](https://blog.codinghorror.com/a-pragmatic-quick-reference/)
you will also find best practices to obtain and maintain good code.

So good code not only just works, it adheres to some basic principles
which we repeat here (using their established buzzwords).

### Use American English

Use it at least for code, comments, documentation and commit messages.

### Code DRY and SOLID

Don't repeat yourself.
Whenever you feel the urge of using copy-and-paste, re-think the design
and use cut-and-paste instead.
*[Don't repeat yourself.](http://c2.com/cgi/wiki?DontRepeatYourself)*

Every class should have a *[single responsibility](https://www.oodesign.com/single-responsibility-principle.html)*, i.e., it does
one thing and only one thing.

Keep your classes *[open for extension but closed for modification](https://www.oodesign.com/open-close-principle.html)*,
i.e., you should be able to write derived classes that extend the behavior of
a class.

Design your classes such that you can always replace an instance
of a base class by any of its derived classes without affecting
correctness
(*[Liskov substitution](https://www.oodesign.com/liskov-s-substitution-principle.html)*).

Whenever you implement a part of an interface by dummy methods (e.g., throwing
NotImplemented exceptions), it is a clear sign that it is not as fine-grained
as it should be.
Apply the principle of *[interfaces segregation](https://www.oodesign.com/interface-segregation-principle.html)*.

High-level classes [should not *depend*](https://www.oodesign.com/dependency-inversion-principle.html)
on low-level classes but on abstractions (interfaces) of low-level classes.

### Write unit tests

A [bandit](https://banditcpp.github.io/bandit/) is waiting for you in the `test` directory.
It wants you to write [unit tests](https://en.wikipedia.org/wiki/Unit_testing).

Note that it may be helpful to write a (failing) test *before* you implement
the actual class. This [development process](https://en.wikipedia.org/wiki/Test-driven_development)
forces you to think about the design, the functionality of the class, corner
cases, and how to test it.

To add a test for a new class, simply write the bandit test code into a new file in the respective subdirectory of `test/src/`.
The file should be named `your-class-name.cpp` or `your-class-category.cpp`.
For test development, it is recommended to set `OGDF_SEPARATE_TESTS=ON` in your CMake configuration.

### Creating new OGDF header or source files

If you create a new header file (extension `.h`),
it belongs into the `include/ogdf/*subdir*` directory.
A new source file (extension `.cpp`) belongs into the
`src/ogdf/*subdir*` directory.
In each case, `*subdir*` is the best-fitting existing subdirectory.

The name of the files (without extension) should match the
name of the main class declared or implemented in that file.
For example, `MyNiftyClass` is declared in `MyNiftyClass.h` and
implemented in `MyNiftyClass.cpp` (unless it is a template).
In any case, start with the [header file](template.h)
and [source file template](template.cpp), respectively.

### Self-documenting code

Comments are very important to describe ideas and explanations,
but best code is self-documenting and easy understandable.
For example, use descriptive identifier (class, method, variable, ...)
names throughout your code.

### "I heard it's faster if I write it in this cryptic way."

This is wrong, most of the time.

If you want to compute `someInt / 2`, your compiler usually knows
that this should be translated to `someInt >> 1`. And it knows
a lot of other tricks. So most of the time it is better to write
what you mean and not the cryptic way that is faster in assembly
language. (If you have a certain exception, move that code into
its own function.)

### Merging branches

We consider 24 hours the minimum timespan between proposing a
merge request and accepting it. If major changes are introduced
to the merge request during this period the review phase will be
extended.

----

## Coding Style

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to ensure
that the code adheres to a certain code style.
Please refer to the [workflow guide](workflow.md#conventions) for more details
on how to use clang-format.
To understand the exact style conventions we use, take a look at the
`.clang-format` file in the main directory; the options are explained by the
[official clang-format documentation](https://clang.llvm.org/docs/ClangFormatStyleOptions.html).
In the following, we give some additional information on style issues that may
not be covered by clang-format.

### On variable names

Prefer camelCase for variable names.
Prefix names of member variables by `m_`, e.g. `m_myValue`.

### On macro names

If you `#define` a macro (even only for temporary use), prefix it by `OGDF_`.
This rule should be applied in general but it must be applied in header files
(with a few basic exceptions) to avoid name clashes with user code and
header files of other libraries.

Write macro names in all-caps with underscore (`_`) as word delimiter.

Also note that double underscores (`__`) must not be used in macro names
since they are reserved for internal use of compilers.

### Do not comment out code

Do not add commented out code.
If you *really* need to, use the `#if 0` ... `#endif`
preprocessor directives (instead of `//` or `/*` ... `*/`).
If you have alternative implementations of some code parts,
we prefer to use preprocessor macros,
i.e., use `#ifdef OGDF_SOME_DESCRIPTIVE_MACRO_NAME` ... `#endif`.
You can then simply enable the code by `#define OGDF_SOME_DESCRIPTIVE_MACRO_NAME`.

### Do not use typedef

Use C++11 `using` instead of `typedef` because it improves readability and
allows templated type definitions.

### Do not use unscoped enumerations

We do not want to use unscoped enumerations.
Those can be safely replaced by either `static const` or
C++11 `enum class`.

----

## Documentation

### Reference documentation

We use [doxygen](https://www.doxygen.org/) to produce the reference documentation
directly from the source code.
Please refer to the doxygen documentation for details about doxygen syntax.

Document (at least) each member in the public interface of the code,
e.g., public methods in a class.

For methods, document the parameters if they are not self-explaining
or not all possible values of a type are a valid input (e.g., if we
have an `unsigned int` but `0` is not allowed).

Each documented member must have a brief comment;
if necessary provide also a long comment and document parameters of functions.

Each class must have at least a brief comment.
Usually classes should have more detailed documentation that also includes
example usage of the class.

We encourage the use of doxygen's special commands whenever referencing
the respective entities. Example:

```c++
/**
 * Returns whether \p e is a bridge.
 *
 * An edge is a bridge iff removing it from the graph increases the number
 * of components.
 *
 * @param G the graph containing \p e.
 * @param e the edge to be tested.
 * @param size is min{ |\a C_u|, |\a C_v| }, where \a C_w is the component
 *        of a node \a w in \p G - \p e, and \p e = (\a u, \a v). This value
 *        is undefined if \p e is not a bridge.
 * @tparam T the type of the returned \p size.
 * @return \c true iff \p e is a bridge of \p G.
 */
template<typename T>
bool isBridge(const Graph &G, edge e, T &size);
```

### Compatibility-breaking changes

If you change the interface of classes that already existed in the last release,
describe these changes in [the Porting Guide](porting/unreleased.md).

### Adding or changing a C++ function macro

If you add or change a function macro (`#define` with parameters) of general interest,
keep [the list of compiler defines and useful macros](defines.md) up-to-date.

### Changes to the build system

If you apply changes to the build system, remember to update the [Build Guide](build.md).

----

## Version Control

We use [git, GitHub, and the respective standard workflow](workflow.md).

### Commit messages

Use sound commit messages.
Use present tense ("Add Foo" not "Added Foo") to describe the changes of
the commit and past tense when you refer to the state before the commit.

The first line of a commit message should be regarded as a subject or title
and is a summary of the commit.
It should not be too long (try to limit it to 80 characters).

The first line can be followed by a blank line and a more detailed description
of the commit, a rationale, etc.

### Commit granularity

Don't do everything in one commit.
You have your own branch where you can have one commit for each
logical unit of changes.

### Naming convention for branches

Use lower-case characters, numbers, and `-` only.
Choose a descriptive name for a feature branch.
Prefix that name with `issue-42-` if the branch is targeted mainly at resolving issue #42.
Then prefix that name with your (or the lead developer's) GitHub username and a `-`.

For example, if Peter wants to get rid of deprecated method calls he could create a branch called `peter-deprecated`.
If deprecated methods were reported as issue #38 beforehand, the branch should instead be called `peter-issue-38-deprecated`.
