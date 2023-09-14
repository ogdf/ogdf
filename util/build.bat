:: A batch script to compile OGDF under windows (for the CI)
::
:: Author: Stephan Beyer

set config=Release
set "cmakeconfig=-DOGDF_SEPARATE_TESTS=1 -DOGDF_WARNING_ERRORS=1"
if "%1" == "debug" (
	set config=Debug
	set "cmakeconfig=%cmakeconfig% -DBUILD_SHARED_LIBS=1"
)

setlocal EnableDelayedExpansion
if "%2" == "gurobi" (
    set "gurobiconfig=-DCOIN_EXTERNAL_SOLVER_LIBRARIES=%GUROBI_HOME%\\lib\\gurobi%GUROBI_VERSION%.lib;%GUROBI_HOME%\\lib\\gurobi_c++"
    if "%1" == "debug" (
        set "gurobiconfig=!gurobiconfig!mdd2019.lib"
    ) else (
        set "gurobiconfig=!gurobiconfig!md2019.lib"
    )

    :: Uses cmake.exe and MSBuild.exe delivered with Visual Studio 2019, found via PATH environment variable
    cmake !gurobiconfig! %cmakeconfig% -DCOIN_SOLVER=GRB -DCOIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES=%GUROBI_HOME%\\include "Visual Studio 16 2019" -A x64 .
) else (
    cmake %cmakeconfig% "Visual Studio 16 2019" -A x64 .
)
endlocal

if %errorlevel% neq 0 (
	exit /b %errorlevel%
)
set msBuildCommand=msbuild OGDF-PROJECT.sln /m /p:Configuration=%config% /nologo /v:m /p:TrackFileAccess=false
%msBuildCommand%
exit /b %errorlevel%
