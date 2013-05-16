# Make VCProj
#
# October 2012
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
	def __init__(self, ttag, tpath, tpats):
		self.tag, self.path, self.pats = ttag, tpath, tpats

def bailout(msg):
	print(msg)
	print('Please use the original makeVCProj.config as a template')
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
			libs += str + ' '
	return libs


#########################################################
# LOAD CONFIGURATION

config = configparser.ConfigParser()

makecfg = 'makeVCProj.config'
for x in sys.argv[1:]:
    if x[:7] == "config=":
        makecfg = x[7:]

print('Loading ' + makecfg + '...')

try:
	config.readfp( open(makecfg))
except IOError:
	bailout(makecfg + ' not found')

if not config.has_section('OGDF'):
	bailout('Section "OGDF" is missing')
if not config.has_section('COIN'):
	bailout('Section "COIN" is missing')

#########################################################
# CONFIGS

# Filenames
filename_vcproj =  loadConfig('OGDF', 'projectFile')
filename_template = loadConfig('OGDF', 'templateFile')
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
	filename_coin_vcproj =  loadConfig('COIN', 'projectFile')
	filename_coin_template = loadConfig('COIN', 'templateFile')

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

	if createDLL and createDLL.startswith('t'):
		addLibs += 'coin.lib '

		addLibsDebugWin32   += getLibs(loadConfig('COIN', 'solverLibs_win32_debug'))
		addLibsReleaseWin32 += getLibs(loadConfig('COIN', 'solverLibs_win32_release'))
		addLibsDebugX64     += getLibs(loadConfig('COIN', 'solverLibs_x64_debug'))
		addLibsReleaseX64   += getLibs(loadConfig('COIN', 'solverLibs_x64_release'))

linkSectionD32 = ''
linkSectionR32 = ''
linkSectionD64 = ''
linkSectionR64 = ''
if createDLL and createDLL.startswith('t'):
	libraryType = '2'
	config_defines += '#define OGDF_DLL\n'
	addDefines += 'OGDF_INSTALL'
	linkerName = 'VCLinkerTool'
	linkSectionBegin = '				AdditionalDependencies=\"psapi.lib '
	linkSectionEnd = '\"\n\
				AdditionalLibraryDirectories="$(SolutionDir)$(PlatformName)\$(ConfigurationName)"\
				LinkTimeCodeGeneration="0"'

	linkSectionD32 = linkSectionBegin + addLibs + addLibsDebugWin32   + linkSectionEnd
	linkSectionR32 = linkSectionBegin + addLibs + addLibsReleaseWin32 + linkSectionEnd
	linkSectionD64 = linkSectionBegin + addLibs + addLibsDebugX64     + linkSectionEnd
	linkSectionR64 = linkSectionBegin + addLibs + addLibsReleaseX64   + linkSectionEnd

else:
	libraryType = '4'
	linkerName = 'VCLibrarianTool'

print('Generating config_autogen.h ...')
config_autogen = open('include\\ogdf\\internal\\config_autogen.h','w')
config_autogen.write('//\n')
config_autogen.write('// This file has been automatically generated by makeVCProj.py\n')
config_autogen.write('//\n')
config_autogen.write('// Do not change this file manually, instead reconfigure OGDF!\n')
config_autogen.write('\n')
config_autogen.write(config_defines)
config_autogen.close()

addIncludes = addIncludes[:-1]

libraryTypeTag = '<<LIBRARYTYPETAG>>'
defineTag = '<<DEFINETAG>>'
includeTag = '<<INCLUDETAG>>'
linkerNameTag = '<<LINKERNAME>>'
linkTagD32 = '<<LINKTAGD32>>'
linkTagR32 = '<<LINKTAGR32>>'
linkTagD64 = '<<LINKTAGD64>>'
linkTagR64 = '<<LINKTAGR64>>'
openMPTag = '<<OPENMPTAG>>'

# Params are:
# - Tag in template-File
# - Directory to start search & subfilters from
# - File Patterns
cppStuff = stuff( '<<CPPTAG>>', 'src\\ogdf', [ '*.c', '*.cpp' ] )
hStuff = stuff( '<<HTAG>>', 'include\\ogdf', [ '*.h' ] )
hLegacyStuff = stuff( '<<HLEGACYTAG>>', 'include\\ogdf_legacy', [ '*.h' ] )
ogdfStuff = [ cppStuff, hStuff, hLegacyStuff ]

cppCoinStuff = stuff( '<<CPPTAG>>', 'src\\coin', [ '*.c', '*.cpp' ] )
hCoinStuff = stuff( '<<HTAG>>', 'include\\coin', [ '*.h', '*.hpp' ] )
coinStuff = [ cppCoinStuff, hCoinStuff ]


#########################################################
#########################################################
## only code below...

# just the def. nothing happens yet
def Walk( curdir, pats, intro ):
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
				if Walk( outpath, pats, intro + '<Filter Name="' + name + '">\n'):
					intro = ''
					vcproj.write('</Filter>\n')
		else:
			for pat in pats:
				if fnmatch.fnmatch(name, pat):
					if len(intro)>0:
						vcproj.write(intro)
						intro = ''
					vcproj.write('<File RelativePath="' + outpath + '"> </File>\n')
	return len(intro) == 0

##########################################
## Main below...

print('Generating VCProj...', end='')

includeLegacyCode = 0;
if len(sys.argv)>1 and sys.argv[1]=='legacy':
	includeLegacyCode = 1
	print('(including legacy code)')

vcproj = open(filename_vcproj,'w')
template = open(filename_template)

check = 0
for line in template:
	if check < len(ogdfStuff) and line.find(ogdfStuff[check].tag) > -1:
		if (ogdfStuff[check].tag!='<<HLEGACYTAG>>' or includeLegacyCode):
			Walk(ogdfStuff[check].path, ogdfStuff[check].pats, '')
		check = check + 1
	elif line.find(libraryTypeTag) > -1:
		vcproj.write(line.replace(libraryTypeTag,libraryType,1))
	elif line.find(openMPTag) > -1:
		vcproj.write(line.replace(openMPTag,ogdfOpenMP,1))
	elif line.find(defineTag) > -1:
		vcproj.write(line.replace(defineTag,addDefines,1))
	elif line.find(includeTag) > -1:
		vcproj.write(line.replace(includeTag,addIncludes,1))
	elif line.find(linkerNameTag) > -1:
		vcproj.write(line.replace(linkerNameTag,linkerName,1))
	elif line.find(linkTagD32) > -1:
		vcproj.write(line.replace(linkTagD32,linkSectionD32,1))
	elif line.find(linkTagR32) > -1:
		vcproj.write(line.replace(linkTagR32,linkSectionR32,1))
	elif line.find(linkTagD64) > -1:
		vcproj.write(line.replace(linkTagD64,linkSectionD64,1))
	elif line.find(linkTagR64) > -1:
		vcproj.write(line.replace(linkTagR64,linkSectionR64,1))
	else:
		vcproj.write(line)

template.close()
vcproj.close()

print('done')

if useCoin:
	print ('Generating COIN VCXProj...', end='')

	vcproj = open(filename_coin_vcproj,'w')
	template = open(filename_coin_template)

	check = 0
	for line in template:
		if check < len(coinStuff) and line.find(coinStuff[check].tag) > -1:
			Walk(coinStuff[check].path, coinStuff[check].pats, '')
			check = check + 1
		elif line.find(openMPTag) > -1:
			vcproj.write(line.replace(openMPTag,coinOpenMP,1))
		elif line.find(includeTag) > -1:
			vcproj.write(line.replace(includeTag,addIncludes,1))
		else:
			vcproj.write(line)

	template.close()
	vcproj.close()

	print ('done')
