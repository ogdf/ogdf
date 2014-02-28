# Make VCXProj
#
# November 2012
# Markus Chimani, markus.chimani@uni-jena.de
# Carsten Gutwenger, carsten.gutwenger@cs.tu-dortmund.de
#########################################################

# use Python 3 print function on both Python 2.x and Python 3.x
from __future__ import print_function

import os, sys, fnmatch

# Simple hack for using "configparser" on both Python 2 and 3
if sys.hexversion > 0x03000000:
	import configparser
else:
	import ConfigParser
	configparser = ConfigParser


class stuff:
	def __init__(self, ttag, tpath, tpats, tcommand, tfilter):
		self.tag, self.path, self.pats, self.command, self.filter = ttag, tpath, tpats, tcommand, tfilter

def bailout(msg):
	print(msg)
	print('Please use the original makeVCXProj.config as a template')
	sys.exit()

def loadConfig(sect, key, noError = False ):
	if config.has_option(sect, key):
		v = config.get(sect, key)
		print('   [', sect, ']', key, '=', v)
		return v
	else:
		if noError:
			return None
		else:
			bailout('Option "' + key + '" in section "' + sect + '" is missing')

def checkSolver(name, defSolver, extSolvers):
	if(name == defSolver):
		return True
	for s in extSolvers:
		s = s.strip()
		if(name == s):
			return True
	return False

def getLibs(libsConfig):
	libs = ''
	for str in libsConfig.split(';'):
		str = str.strip();
		if str != '':
			libs += str + ';'
	return libs


#########################################################
# LOAD CONFIGURATION

config = configparser.ConfigParser()

makecfg = 'makeVCXProj.config'
for x in sys.argv[1:]:
	if x[:7] == "config=":
		makecfg = x[7:]

print('Loading ' + makecfg + '...')

try:
	config.readfp( open(makecfg))
except IOError:
	bailout(makecfg + ' not found')

if not config.has_section('GENERAL'):
	bailout('Section "GENERAL" is missing')
if not config.has_section('OGDF'):
	bailout('Section "OGDF" is missing')
if not config.has_section('OGDF-TEST'):
	bailout('Section "OGDF-TEST" is missing')
if not config.has_section('COIN'):
	bailout('Section "COIN" is missing')

#########################################################
# CONFIGS

platformToolset = loadConfig('GENERAL', 'platformToolset')
windowsVersion = loadConfig('GENERAL', 'windowsVersion')
createSolution = loadConfig('GENERAL', 'createSolution').startswith('t')

# Filenames
filename_sln = loadConfig('GENERAL', 'solutionFile')
filename_vcxproj =  loadConfig('OGDF', 'projectFile')
filename_template = loadConfig('OGDF', 'templateFile')
filename_vcxfilters =  loadConfig('OGDF', 'projectFiltersFile')
filename_template_filters = loadConfig('OGDF', 'templateFiltersFile')
addDefines = ''
addIncludes = ''
addLibs = ''
addLibsWin32 = ''
addLibsWin32 = ''
addLibsX64 = ''
addLibsX64 = ''

# defines for auto generated config_autogen.h header file
config_defines = ''

createDLL = loadConfig('OGDF', 'DLL', True)

memoryManager = loadConfig('OGDF', 'memoryManager', True)
if memoryManager:
	config_defines += '#define ' + memoryManager + '\n'

ogdfOpenMP = 'false'
coinOpenMP = 'false'
ogdfEnableOpenMP = loadConfig('OGDF', 'OpenMP', True)
if ogdfEnableOpenMP and ogdfEnableOpenMP.startswith('t'):
	ogdfOpenMP = 'true'

useCoin = loadConfig('COIN', 'useCoin').startswith('t')
addOsiCpx = False
addOsiGrb = False
addLibsDebugWin32 = ''
addLibsReleaseWin32 = ''
addLibsDebugX64 = ''
addLibsReleaseX64 = ''
if useCoin:
	filename_coin_vcxproj =  loadConfig('COIN', 'projectFile')
	filename_coin_template = loadConfig('COIN', 'templateFile')
	filename_coin_vcxfilters =  loadConfig('COIN', 'projectFiltersFile')
	filename_coin_template_filters = loadConfig('COIN', 'templateFiltersFile')

	coinEnableOpenMP = loadConfig('COIN', 'OpenMP', True)
	if coinEnableOpenMP and coinEnableOpenMP.startswith('t'):
		coinOpenMP = 'true'

	defaultSolver   = loadConfig('COIN', 'defaultSolver')
	externalSolvers = loadConfig('COIN', 'externalSolvers').split(';')
	solverIncludes  = loadConfig('COIN', 'solverIncludes').split(';')

	addOsiCpx = checkSolver('CPX', defaultSolver, externalSolvers)
	addOsiGrb = checkSolver('GRB', defaultSolver, externalSolvers)

	config_defines += '#define USE_COIN\n'
	config_defines += '#define COIN_OSI_' + defaultSolver + '\n'

	for s in externalSolvers:
		s = s.strip();
		if s != '':
			config_defines += '#define OSI_' + s.strip() + '\n'

	for p in solverIncludes:
		p = p.strip()
		if p != '':
			addIncludes += p + ';'

	addLibs += 'coin.lib;'

	addLibsDebugWin32   += getLibs(loadConfig('COIN', 'solverLibs_win32_debug'))
	addLibsReleaseWin32 += getLibs(loadConfig('COIN', 'solverLibs_win32_release'))
	addLibsDebugX64     += getLibs(loadConfig('COIN', 'solverLibs_x64_debug'))
	addLibsReleaseX64   += getLibs(loadConfig('COIN', 'solverLibs_x64_release'))

linkSectionD32 = ''
linkSectionR32 = ''
linkSectionD64 = ''
linkSectionR64 = ''
if createDLL and createDLL.startswith('t'):
	libraryType = 'DynamicLibrary'
	config_defines += '#define OGDF_DLL\n'
	addDefines += 'OGDF_INSTALL'

	linkSectionBegin = '    <Link>\n\
	<AdditionalDependencies>';
	linkSectionEnd = 'kernel32.lib;user32.lib;gdi32.lib;winspool.lib;comdlg32.lib;advapi32.lib;shell32.lib;psapi.lib;%(AdditionalDependencies)</AdditionalDependencies>\n\
	<LinkTimeCodeGeneration>\n\
	</LinkTimeCodeGeneration>\n\
    <AdditionalLibraryDirectories>$(SolutionDir)$(Platform)\$(Configuration)\</AdditionalLibraryDirectories>\n\
    </Link>'

	linkSectionD32 = linkSectionBegin + addLibs + addLibsDebugWin32   + linkSectionEnd
	linkSectionR32 = linkSectionBegin + addLibs + addLibsReleaseWin32 + linkSectionEnd
	linkSectionD64 = linkSectionBegin + addLibs + addLibsDebugX64     + linkSectionEnd
	linkSectionR64 = linkSectionBegin + addLibs + addLibsReleaseX64   + linkSectionEnd

else:
	libraryType = 'StaticLibrary'

linkSectionProgramBegin = '    <Link>\n\
    <SubSystem>Console</SubSystem>\n\
    <AdditionalDependencies>';
linkSectionProgramEnd = 'ogdf.lib;kernel32.lib;user32.lib;gdi32.lib;winspool.lib;comdlg32.lib;advapi32.lib;shell32.lib;psapi.lib;%(AdditionalDependencies)</AdditionalDependencies>\n\
    <AdditionalLibraryDirectories>$(SolutionDir)$(Platform)\$(Configuration)\</AdditionalLibraryDirectories>\n\
    </Link>'

linkSectionProgramD32 = linkSectionProgramBegin + addLibs + addLibsDebugWin32   + linkSectionProgramEnd
linkSectionProgramR32 = linkSectionProgramBegin + addLibs + addLibsReleaseWin32 + linkSectionProgramEnd
linkSectionProgramD64 = linkSectionProgramBegin + addLibs + addLibsDebugX64     + linkSectionProgramEnd
linkSectionProgramR64 = linkSectionProgramBegin + addLibs + addLibsReleaseX64   + linkSectionProgramEnd

	
filename_ogdf_test_vcxproj =  loadConfig('OGDF-TEST', 'projectFile')
filename_ogdf_test_template = loadConfig('OGDF-TEST', 'templateFile')
filename_ogdf_test_vcxfilters =  loadConfig('OGDF-TEST', 'projectFiltersFile')
filename_ogdf_test_template_filters = loadConfig('OGDF-TEST', 'templateFiltersFile')

print('Generating config_autogen.h ...')
config_autogen = open('include\\ogdf\\internal\\config_autogen.h','w')
config_autogen.write('//\n')
config_autogen.write('// This file has been automatically generated by makeVCXProj.py\n')
config_autogen.write('//\n')
config_autogen.write('// Do not change this file manually, instead reconfigure OGDF!\n')
config_autogen.write('\n')
config_autogen.write(config_defines)
config_autogen.write('\n#ifndef _WIN32_WINNT\n')
config_autogen.write('#define _WIN32_WINNT ' + windowsVersion + '\n')
config_autogen.write('#endif\n')
config_autogen.close()

addIncludes = addIncludes[:-1]
libraryTypeTag = '<<LIBRARYTYPETAG>>'
defineTag = '<<DEFINETAG>>'
includeTag = '<<INCLUDETAG>>'
linkTagD32 = '<<LINKTAGD32>>'
linkTagR32 = '<<LINKTAGR32>>'
linkTagD64 = '<<LINKTAGD64>>'
linkTagR64 = '<<LINKTAGR64>>'
filtersTag = '<<FTAG>>'
openMPTag = '<<OPENMPTAG>>'
toolsetTag = '<<TOOLSET>>'

# Params are:
# - Tag in template-File
# - Directory to start search & subfilters from
# - File Patterns
cppStuff = stuff( '<<CPPTAG>>', 'src\\ogdf', [ '*.c', '*.cpp' ], 'ClCompile', 'Source Files' )
hStuff = stuff( '<<HTAG>>', 'include\\ogdf', [ '*.h', '*.inc' ], 'ClInclude', 'Header Files' )
hLegacyStuff = stuff( '<<HLEGACYTAG>>', 'include\\ogdf_legacy', [ '*.h' ], 'ClInclude', 'Header Files Legacy' )
ogdfStuff = [ cppStuff, hStuff, hLegacyStuff ]

cppCoinStuff = stuff( '<<CPPTAG>>', 'src\\coin', [ '*.c', '*.cpp' ], 'ClCompile', 'Source Files' )
hCoinStuff = stuff( '<<HTAG>>', 'include\\coin', [ '*.h', '*.hpp' ], 'ClInclude', 'Header Files' )
coinStuff = [ cppCoinStuff, hCoinStuff ]

cppTestStuff = stuff( '<<CPPTAG>>', 'test', [ '*.c', '*.cpp' ], 'ClCompile', 'Source Files' )
hTestStuff = stuff( '<<HTAG>>', 'test', [ '*.h', '*.hpp' ], 'ClInclude', 'Header Files' )
testStuff = [ cppTestStuff, hTestStuff ]

#########################################################
#########################################################
## only code below...

# just the def. nothing happens yet
def Walk( curdir, pats, command ):
	names = os.listdir( curdir)
	names.sort()

	for name in names:
		# OGDF ignores
		if name.startswith('.') or name.startswith('_') or (name=='legacy' and not includeLegacyCode):
			continue
		# COIN ignores
		if (name == 'OsiCpxSolverInterface.cpp' or name == 'OsiCpxSolverInterface.hpp') and not addOsiCpx:
			continue
		if (name == 'OsiGrbSolverInterface.cpp' or name == 'OsiGrbSolverInterface.hpp') and not addOsiGrb:
			continue

		outpath = curdir + '\\' + name
		fullname = os.path.normpath(outpath)

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			if name != 'abacus' or useCoin == True:
				Walk( outpath, pats, command)
		else:
			for pat in pats:
				if fnmatch.fnmatch(name, pat):
					vcxproj.write('    <' + command + ' Include="' + outpath + '" />\n')

def WalkFilterFiles( curdir, pats, command, filter ):
	names = os.listdir( curdir)
	names.sort()

	for name in names:
		# OGDF ignores
		if name.startswith('.') or name.startswith('_') or (name=='legacy' and not includeLegacyCode):
			continue
		# COIN ignores
		if (name == 'OsiCpxSolverInterface.cpp' or name == 'OsiCpxSolverInterface.h') and not addOsiCpx:
			continue
		if (name == 'OsiGrbSolverInterface.cpp' or name == 'OsiGrbSolverInterface.h') and not addOsiGrb:
			continue

		outpath = curdir + '\\' + name
		fullname = os.path.normpath(outpath)

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			WalkFilterFiles( outpath, pats, command, filter + '\\' + name )
		else:
			for pat in pats:
				if fnmatch.fnmatch(name, pat):
					vcxfilters.write('    <' + command + ' Include="' + outpath + '">\n')
					vcxfilters.write('      <Filter>' + filter + '</Filter>\n')
					vcxfilters.write('    </' + command + '>\n')

def WalkFilters( curdir, filter ):
	names = os.listdir( curdir)
	names.sort()

	for name in names:
		if name.startswith('.') or name.startswith('_') or (name=='legacy' and not includeLegacyCode):
			continue

		outpath = curdir + '\\' + name
		fullname = os.path.normpath(outpath)

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			if name != 'abacus' or useCoin == True:
				filtername = filter + '\\' + name
				vcxfilters.write('    <Filter Include="' + filtername + '">\n')
				vcxfilters.write('    </Filter>\n')
				WalkFilters( outpath, filtername )


##########################################
## Main below...

print('Generating OGDF VCXProj...', end='')

includeLegacyCode = 0;
if len(sys.argv)>1 and sys.argv[1]=='legacy':
	includeLegacyCode = 1
	print('(including legacy code)...'),

vcxproj = open(filename_vcxproj,'w')
template = open(filename_template)

check = 0
for line in template:
	if check < len(ogdfStuff) and line.find(ogdfStuff[check].tag) > -1:
		if (ogdfStuff[check].tag!='<<HLEGACYTAG>>' or includeLegacyCode):
			Walk(ogdfStuff[check].path, ogdfStuff[check].pats, ogdfStuff[check].command)
		check = check + 1
	elif line.find(libraryTypeTag) > -1:
		vcxproj.write(line.replace(libraryTypeTag,libraryType,1))
	elif line.find(openMPTag) > -1:
		vcxproj.write(line.replace(openMPTag,ogdfOpenMP,1))
	elif line.find(defineTag) > -1:
		vcxproj.write(line.replace(defineTag,addDefines,1))
	elif line.find(includeTag) > -1:
		vcxproj.write(line.replace(includeTag,addIncludes,1))
	elif line.find(linkTagD32) > -1:
		vcxproj.write(line.replace(linkTagD32,linkSectionD32,1))
	elif line.find(linkTagR32) > -1:
		vcxproj.write(line.replace(linkTagR32,linkSectionR32,1))
	elif line.find(linkTagD64) > -1:
		vcxproj.write(line.replace(linkTagD64,linkSectionD64,1))
	elif line.find(linkTagR64) > -1:
		vcxproj.write(line.replace(linkTagR64,linkSectionR64,1))
	elif line.find(toolsetTag) > -1:
		vcxproj.write(line.replace(toolsetTag,platformToolset,1))
	else:
		vcxproj.write(line)

template.close()
vcxproj.close()

# Creation of filters file...

vcxfilters = open(filename_vcxfilters,'w')
template_filters = open(filename_template_filters)

check = 0
for line in template_filters:
	if check < len(ogdfStuff) and line.find(ogdfStuff[check].tag) > -1:
		if (ogdfStuff[check].tag!='<<HLEGACYTAG>>' or includeLegacyCode):
			WalkFilterFiles(ogdfStuff[check].path, ogdfStuff[check].pats, ogdfStuff[check].command, ogdfStuff[check].filter)
		check = check + 1
	elif line.find(filtersTag) > -1:
		for s in ogdfStuff:
			if (s.tag!='<<HLEGACYTAG>>' or includeLegacyCode):
				WalkFilters(s.path, s.filter)
	else:
		vcxfilters.write(line)

template_filters.close()
vcxfilters.close()

print('done')

if useCoin:
	print ('Generating COIN VCXProj...', end='')

	vcxproj = open(filename_coin_vcxproj,'w')
	template = open(filename_coin_template)

	check = 0
	for line in template:
		if check < len(coinStuff) and line.find(coinStuff[check].tag) > -1:
			Walk(coinStuff[check].path, coinStuff[check].pats, coinStuff[check].command)
			check = check + 1
		elif line.find(openMPTag) > -1:
			vcxproj.write(line.replace(openMPTag,coinOpenMP,1))
		elif line.find(includeTag) > -1:
			vcxproj.write(line.replace(includeTag,addIncludes,1))
		elif line.find(toolsetTag) > -1:
			vcxproj.write(line.replace(toolsetTag,platformToolset,1))
		else:
			vcxproj.write(line)

	template.close()
	vcxproj.close()

	# Creation of filters file...

	vcxfilters = open(filename_coin_vcxfilters,'w')
	template_filters = open(filename_coin_template_filters)

	check = 0
	for line in template_filters:
		if check < len(coinStuff) and line.find(coinStuff[check].tag) > -1:
			WalkFilterFiles(coinStuff[check].path, coinStuff[check].pats, coinStuff[check].command, coinStuff[check].filter)
			check = check + 1
		elif line.find(filtersTag) > -1:
			for s in coinStuff:
				WalkFilters(s.path, s.filter)
		else:
			vcxfilters.write(line)

	template_filters.close()
	vcxfilters.close()

	print ('done')

print ('Generating OGDF-TEST VCXProj...', end='')

vcxproj = open(filename_ogdf_test_vcxproj,'w')
template = open(filename_ogdf_test_template)

check = 0
for line in template:
	if check < len(testStuff) and line.find(testStuff[check].tag) > -1:
		Walk(testStuff[check].path, testStuff[check].pats, testStuff[check].command)
		check = check + 1
	elif line.find(includeTag) > -1:
		vcxproj.write(line.replace(includeTag,addIncludes,1))
	elif line.find(linkTagD32) > -1:
		vcxproj.write(line.replace(linkTagD32,linkSectionProgramD32,1))
	elif line.find(linkTagR32) > -1:
		vcxproj.write(line.replace(linkTagR32,linkSectionProgramR32,1))
	elif line.find(linkTagD64) > -1:
		vcxproj.write(line.replace(linkTagD64,linkSectionProgramD64,1))
	elif line.find(linkTagR64) > -1:
		vcxproj.write(line.replace(linkTagR64,linkSectionProgramR64,1))
	elif line.find(toolsetTag) > -1:
		vcxproj.write(line.replace(toolsetTag,platformToolset,1))
	else:
		vcxproj.write(line)

template.close()
vcxproj.close()

# Creation of filters file...

vcxfilters = open(filename_ogdf_test_vcxfilters,'w')
template_filters = open(filename_ogdf_test_template_filters)

check = 0
for line in template_filters:
	if check < len(testStuff) and line.find(testStuff[check].tag) > -1:
		WalkFilterFiles(testStuff[check].path, testStuff[check].pats, testStuff[check].command, testStuff[check].filter)
		check = check + 1
	elif line.find(filtersTag) > -1:
		for s in testStuff:
			WalkFilters(s.path, s.filter)
	else:
		vcxfilters.write(line)

template_filters.close()
vcxfilters.close()

print ('done')

if createSolution:
	print ('Generating Solution...', end='')

	GUID_sln = '{0F7C385F-D08C-494E-8715-70CD736C75B2}'
	GUID_ogdf = '{7801D1BE-E2FE-476B-A4B4-5D27F387F479}'
	GUID_coin = '{FB212DCC-D374-430B-B594-5CEC25BC7A75}'
	GUID_ogdf_test = '{278A54D6-588A-4110-9E44-DA24DB9D5E85}'

	Configurations = [ 'Debug', 'Release' ]
	Platforms = [ 'Win32', 'x64' ]

	sln = open(filename_sln, 'w')
	sln.write('Microsoft Visual Studio Solution File, Format Version 11.00\n')
	sln.write('Project("' + GUID_sln + '") = "ogdf", "' + filename_vcxproj + '", "' + GUID_ogdf + '"\n')
	if createDLL and createDLL.startswith('t') and useCoin:
		sln.write('\tProjectSection(ProjectDependencies) = postProject\n')
		sln.write('\t\t' + GUID_coin + ' = ' + GUID_coin + '\n')
		sln.write('\tEndProjectSection\n')
	sln.write('EndProject\n')

	if useCoin:
		sln.write('Project("' + GUID_sln + '") = "coin", "' + filename_coin_vcxproj + '", "' + GUID_coin + '"\n')
		sln.write('EndProject\n')

	sln.write('Project("' + GUID_sln + '") = "ogdf-test", "' + filename_ogdf_test_vcxproj + '", "' + GUID_ogdf_test + '"\n')
	sln.write('EndProject\n')

	sln.write('Global\n')
	sln.write('\tGlobalSection(SolutionConfigurationPlatforms) = preSolution\n')
	sln.write('\t\tDebug|Win32 = Debug|Win32\n')
	sln.write('\t\tDebug|x64 = Debug|x64\n')
	sln.write('\t\tRelease|Win32 = Release|Win32\n')
	sln.write('\t\tRelease|x64 = Release|x64\n')
	sln.write('\tEndGlobalSection\n')
	sln.write('\tGlobalSection(ProjectConfigurationPlatforms) = postSolution\n')

	if useCoin:
		projectGUIDs = [ GUID_ogdf, GUID_coin, GUID_ogdf_test ]
	else:
		projectGUIDs = [ GUID_ogdf, GUID_ogdf_test ]
	for g in projectGUIDs:
		for c in Configurations:
			for p in Platforms:
				sln.write('\t\t' + g + '.' + c + '|' + p + '.ActiveCfg' + ' = ' + c + '|' + p + '\n')
				sln.write('\t\t' + g + '.' + c + '|' + p + '.Build.0' + ' = ' + c + '|' + p + '\n')
	sln.write('\tEndGlobalSection\n')
	sln.write('\tGlobalSection(SolutionProperties) = preSolution\n')
	sln.write('\t\tHideSolutionNode = FALSE\n')
	sln.write('\tEndGlobalSection\n')
	sln.write('EndGlobal\n')

	sln.close()
	print ('done')
