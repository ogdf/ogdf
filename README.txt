***********************************************************
*         OGDF - The Open Graph Drawing Framework         *
*                                                         *
*                         README                          *
***********************************************************

Welcome to OGDF!

OGDF is a portable C++ class library for graph drawing.
This archive contains the source-code of OGDF.


******************** LICENSE ********************

This software is distributed under the terms of the GNU
General Public License v2 or v3, with special exceptions.
By installing this software you agree to these license terms.

See LICENSE.txt for more information on the license and
included third-party software, like frameworks for
linear programming or unit testing.

If you have questions, please contact us.


******************* COPYRIGHT *******************

All files in the OGDF distribution (except for third-party
software) are copyrighted:

Copyright (C) 2005-2015


****************** REQUIREMENTS *****************

To be able to compile OGDF, you need a C++11 compiler.

The following compilers should work:
  * Microsoft Visual C++ 2015
  * LLVM/clang >= 3.5
  * G++ >= 4.8

It is known that Visual C++ 2014, G++ 4.7, LLVM/clang 3.2
do not support a sufficient subset of the C++11 standard
to compile OGDF.

To generate the build files for OGDF, you can choose
between the "old way" or using CMake >= 3.1.
For the "old way", you need Python >= 2.6.
We are currently migrating to CMake, so problems
might occur.

You need Doxygen >= 1.8 to build the documentation.


****************** INSTALLATION *****************


Build OGDF using CMake [any platform]:

  CMake can generate Makefiles as well as project files
  for many IDEs (including Visual Studio and Xcode).
  You can use CMake the standard way. Running cmake with
  the source path and some additional options (like the
  generator to use) should be sufficient, then proceed
  as usual with your IDE or `make`.

  There are also GUIs (cmake-gui) and UIs (ccmake) available
  where you can easily configure (and generate) the build
  system.
  Note that we have defined OGDF-specific options, for example,
  to specify an external LP solver to be used.

Please refer to the OGDF Wiki for more detailed information:

http://www.ogdf.net/ogdf.php/tech:build


Build OGDF the old way [Linux, Mac OS]:

  1. Edit makeMakefile.config for your configuration
     (if necessary): check the [GENERAL] section. If
     you do not want to use CPLEX or Gurobi LP solvers,
     the default parameters should be suitable.

  2. Execute makeMakefile.py to generate a suitable Makefile.

  3. Call `make` or `make release` to build the OGDF
     library. The compiled and linked files can be found
     in the _release directory.
     You can also run `make debug` to build a debuggable
     version of the OGDF library (including a lot of
     further checks, so it is much slower), whose output
     files can be found in the _debug directory.

  4. Call `make doc` to generate the documentation.


Build OGDF the old way [Microsoft Visual C++ 2015]:

  1. Edit makeVCXProj.config for your configuration
     (if necessary).

  2. Execute the Python script makeVCXProj.py to
     generate the Visual Studio 2015 project file.

  3. Open the created project file ogdf.vcxproj with
     Visual Studio and call build.


******************** CHANGES ********************

For changes refer to the version history at

http://www.ogdf.net/doku.php?id=tech:versions


******************** CONTACT ********************

Web:          http://www.ogdf.net
Mailing list: https://groups.google.com/group/ogdf


Enjoy!

  The OGDF Team.
