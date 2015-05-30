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


****************** INSTALLATION *****************

Unpack the OGDF archive in the directory, where you want to
install OGDF.

Build OGDF [Linux, Mac OS]:

  1. Edit makeMakefile.config for your configuration
     (if necessary): check the [GENERAL] section. If
     you do not want to use CPLEX or Gurobi LP solvers,
     the default parameters should be suitable.

  2. Execute makeMakefile.sh to generate a suitable Makefile.

     MAC USERS: If step 2 fails, go back to step 1 and use
     the alternative compilerCommand configuration!

  3. Call make to build the OGDF library. You may also call
     make debug to generate a debuggable version.


Build OGDF [Microsoft Visual Studio]:

  1. Create Visual Studio project file:

     Visual Studio 2008: Execute the python script makeVCProj.py
       to generate a Visual Studio 2008 project file ogdf.vcproj.

     Visual Studio 2010 (2012, 2013):
           Execute the python script makeVCXProj.py to generate a
           Visual Studio 2010 project file ogdf.vcxproj.
           For VS 2012, 2013, you can edit makeVCXProj.config and
           set the compiler version.

  2. Open the created project file (.vcproj or .vcxproj) with
     Visual Studio and call build.


Please refer to the OGDF Wiki for more detailed information:

gcc:        http://www.ogdf.net/ogdf.php/tech:installgcc
Visual C++: http://www.ogdf.net/ogdf.php/tech:installvcc


******************** CHANGES ********************

For changes refer to the version history at

http://www.ogdf.net/doku.php?id=tech:versions


******************** CONTACT ********************

Web:          http://www.ogdf.net
Mailing list: https://groups.google.com/group/ogdf


Enjoy!

  The OGDF Team.
