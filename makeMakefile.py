#!/usr/bin/env python
# Make Makefile
#
# October 2012
# Markus Chimani, markus.chimani@cs.tu-dortmund.de
# Carsten Gutwenger, carsten.gutwenger@cs.tu-dortmund.de
#########################################################

# use Python 3 print function on both Python 2.x and Python 3.x
from __future__ import print_function

import os, sys, fnmatch, posixpath, shutil

# Simple hack for using "configparser" on both Python 2 and 3
if sys.hexversion > 0x03000000:
	import configparser
else:
	import ConfigParser
	configparser = ConfigParser


class versionclass:
	def call(self):
		return '$(' + self.var + ')'
	def library(self):
		return self.call() + '/' + libName
	def sharedlibrary(self):
		return self.call() + '/' + sharedlibName
	def objects(self):
		return '$(' +self.var + '_OBJS)'
	def path(self):
		return '_' + self.var
	def coinLibrary(self):
		return self.call() + '/' + coinLibName
	def coinSharedLibrary(self):
		return self.call() + '/' + coinSharedLibName
	def coinObjects(self):
		return '$(' +self.var + '_COIN_OBJS)'

def bailout(msg):
	print(msg)
	print('Please use the original makeMakefile.config as a template')
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


#########################################################
# LOAD CONFIGURATION

config = configparser.ConfigParser()
print('Loading makeMakefile.config...')

try:
	config.readfp( open('makeMakefile.config') )
except IOError:
	bailout('makeMakefile.config not found')

if not config.has_section('GENERAL'):
	bailout('Section "GENERAL" is missing')
if not config.has_section('OGDF'):
	bailout('Section "GENERAL" is missing')
if not config.has_section('VERSIONS'):
	bailout('Section "VERSIONS" is missing')
if not config.has_section('COIN'):
	bailout('Section "COIN" is missing')

compilerCommand = loadConfig('GENERAL', 'compilerCommand')
compilerParams = loadConfig('GENERAL', 'compilerParams')
libCommand = loadConfig('GENERAL', 'libCommand')
sharedlibCommand = loadConfig('GENERAL', 'sharedlibCommand')
rmCommand = loadConfig('GENERAL', 'rmCommand')
mkdirCommand = loadConfig('GENERAL', 'mkdirCommand')
ranlibCommand = loadConfig('GENERAL', 'ranlibCommand', True)
if ranlibCommand == None:
	ranlibCommand = ''
installPrefix = loadConfig('GENERAL', 'installPrefix', True)

sharedLib = loadConfig('GENERAL', 'sharedLib').startswith('t')
libName = loadConfig('OGDF', 'libName')
sharedlibName = loadConfig('OGDF', 'sharedlibName')
includeLegacyCode = loadConfig('OGDF', 'includeLegacyCode').startswith('t')
memoryManager = loadConfig('OGDF', 'memoryManager', True)

# defines for auto generated config_autogen.h header file
config_defines = ''

gccMessageLength = loadConfig('GENERAL', 'gccMessageLength', True)
if gccMessageLength == None:
	gccMessageLength = ''
else:
	gccMessageLength = '-fmessage-length=' + gccMessageLength

compiler = ' '.join( [ compilerCommand, gccMessageLength, compilerParams ] )

libs = ''

if sharedLib:
	config_defines += '#define OGDF_DLL\n'
	if sys.platform == 'win32' or sys.platform == 'cygwin':
		libs = ' '.join( [libs, '-lpsapi'] )
	else:
		compiler = ' '.join( [compiler, '-fPIC'] )

if memoryManager:
	config_defines += '#define ' + memoryManager + '\n'

addIncludes = ''
useCoin = loadConfig('COIN', 'useCoin').startswith('t')
addOsiCpx = False
addOsiGrb = False
if useCoin:
	shutil.copyfile('config/coinstuff/config.h', 'include/coin/config.h')
	coinLibName = loadConfig('COIN', 'libName')
	coinSharedLibName = loadConfig('COIN', 'sharedlibName')

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
			addIncludes += '-I' + p + ' '

	if sharedLib and (sys.platform == 'win32' or sys.platform == 'cygwin'):
		libs = ' '.join( [libs, '-lCOIN'] )

ogdfCompiler = compiler + ' -I./include ' + addIncludes
coinCompiler = compiler + ' $(COIN_INSTALL_DEFINES) -I./include/coin ' + addIncludes

if sharedLib:
	ogdfCompiler = ' '.join( [ogdfCompiler, '-DOGDF_INSTALL ' ] )


versions = []
V = config.items('VERSIONS')
if len(V) == 0:
	bailout('Versions missing')
else:
	for ve in V:
		v = versionclass()
		v.var, v.params = ve
		print('   [ VERSIONS ] Name:', v.var, ', Cmd:',v.params)
		versions.append(v)

print('Resulting compiler call (OGDF):', ogdfCompiler)
if useCoin:
	print('Resulting compiler call (COIN):', coinCompiler)

print('Finished loading makeMakefile.config')

#########################################################
# ANALYZE & GENERATE

print('Generating config_autogen.h ...')
config_autogen = open('include/ogdf/internal/config_autogen.h','w')
config_autogen.write('//\n')
config_autogen.write('// This file has been automatically generated by makeMakefile.py\n')
config_autogen.write('//\n')
config_autogen.write('// Do not change this file manually, instead reconfigure OGDF!\n')
config_autogen.write('\n')
config_autogen.write(config_defines)
config_autogen.close()

if installPrefix:
	print('Generating ogdf.pc ...')
	pc = open('ogdf.pc', 'w')
	pc.write('prefix=' + installPrefix + '\n')
	pc.write('exec_prefix=${prefix}\n')
	pc.write('includedir=${prefix}/include\n')
	pc.write('libdir=${exec_prefix}/lib\n\n')
	pc.write('Name: OGDF\n')
	pc.write('Description: Open Graph Drawing Framework\n')
	pc.write('URL: http://ogdf.net/\n')
	pc.write('Version: ')
	version_h = open('include/ogdf/internal/version.h')
	for line in version_h:
		if line.startswith('#define OGDF_VERSION '):
			pc.write(line[22:line.find('"', 22)])
	version_h.close()
	pc.write('\n')
	pc.write('Cflags: -I${includedir}\n')
	pc.write('Libs: -L${libdir} -lOGDF')
	if useCoin:
		pc.write(' -lCOIN')
	pc.write('\n')
	pc.close()

print('Analyzing sources & generating Makefile...')

makefile = open('Makefile','w')

# add header
header = open('Makefile.header')
headercontent = header.read()
header.close()
makefile.write(headercontent)

# define release & debug

for v in versions:
	makefile.write(v.var + ' = ' + v.path() + '\n')
makefile.write('\n');

# just the def. nothing happens yet
def Walk( curdir ):

	objs = []
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

		fullname = posixpath.normpath(posixpath.join(curdir, name))

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			if name != 'abacus' or useCoin == True:
				objs = objs + Walk( fullname )
		else:
			for pat in [ '*.c', '*.cpp' ]:
				if fnmatch.fnmatch(name, pat):
					objfullname = fullname[:-len(pat)+2] + 'o'
					objs.append(objfullname)

					callForDeps = ogdfCompiler + '-MM ' + fullname + ' > targetAndDepend'
					os.system( callForDeps )
					t = open('targetAndDepend')
					targetAndDepend = t.read()
					t.close()

					for v in versions:
						# print target&depend: add full path spec, incl. version & ignore extra line
						path = v.call() + '/' +fullname[:-len(name)]
						makefile.write(path + targetAndDepend[:-1] + '\n')

						# ensure folder
						makefile.write('\t' + mkdirCommand + ' ' + v.call() + '/' + fullname[:-len(name)-1] + '\n')
						# what to do: call the compiler
						makefile.write('\t' + ogdfCompiler + ' ' + v.params + ' -o ' + v.call() + '/' + objfullname + ' -c ' + fullname + '\n\n')

					# pattern found: don't try other suffix
					break
	return objs

def WalkCoin( curdir ):

	objs = []
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

		fullname = posixpath.normpath(posixpath.join(curdir, name))

		if os.path.isdir(fullname) and not os.path.islink(fullname):
			if name != 'abacus' or useCoin == True:
				objs = objs + WalkCoin( fullname )
		else:
			for pat in [ '*.c', '*.cpp' ]:
				if fnmatch.fnmatch(name, pat):
					objfullname = fullname[:-len(pat)+2] + 'o'
					if not name.startswith('unitTest'): #conflict if in library!
						objs.append(objfullname)

					#callForDeps = coinCompiler + '-MM ' + fullname + ' > targetAndDepend'
					#print callForDeps
					#os.system( callForDeps )
					#t = open('targetAndDepend')
					#targetAndDepend = t.read()
					#t.close()

					for v in versions:
						# print target&depend: add full path spec, incl. version & ignore extra line
						#path = v.call() + '/' + fullname[:-len(name)]
						path = v.call() + '/' +objfullname
						#makefile.write(path + targetAndDepend[:-1] + '\n')
						makefile.write(path + ': ' + fullname + '\n')

						# ensure folder
						makefile.write('\t' + mkdirCommand + ' ' + v.call() + '/' + fullname[:-len(name)-1] + '\n')
						# what to do: call the compiler
						makefile.write('\t' + coinCompiler + ' ' + v.params + ' -o ' + v.call() + '/' + objfullname + ' -c ' + fullname + '\n\n')

					# pattern found: don't try other suffix
					break
	return objs


# OGDF
objs = Walk( './src/ogdf' )
# Clean up
os.system(rmCommand + ' targetAndDepend')

# List all Objs for use in lib-generation and clear
for v in versions:
	makefile.write(v.objects()[2:-1] + ' = \\\n')
	for o in objs:
		makefile.write(v.call() + '/' + o + ' \\\n')
	makefile.write('\n')

# Coin
if useCoin:
	objsCoin = WalkCoin( './src/coin' )
	# Clean up
	#os.system(rmCommand + ' targetAndDepend')

	# List all Objs for use in lib-generation and clear
	for v in versions:
		makefile.write(v.coinObjects()[2:-1] + ' = \\\n')
		for o in objsCoin:
			makefile.write(v.call() + '/' + o + ' \\\n')
		makefile.write('\n')

# generate alls and cleans etc...

for v in versions:

	makefile.write('\n#######################################################')
	makefile.write('\n# all, clean, etc. for ' + v.var + '\n\n')

	if sharedLib:
		if useCoin:
			makefile.write(v.coinSharedLibrary() + ': ' + v.coinObjects() + '\n')
			makefile.write('\t' + sharedlibCommand  + ' -shared -o ' + v.coinSharedLibrary() + ' ' + v.coinObjects() + '\n\n')

		makefile.write(v.sharedlibrary() + ': ' + v.objects() + '\n')
		makefile.write('\t' + sharedlibCommand  + ' -shared -o ' + v.sharedlibrary() + ' ' + v.objects() + ' ' + libs + ' $(LIBS)\n\n')

		makefile.write(v.var + ': ' + v.sharedlibrary())
		if useCoin:
			makefile.write(' ' + v.coinSharedLibrary())
		makefile.write('\n\n')

	else:
		if useCoin:
			makefile.write(v.coinLibrary() + ': ' + v.coinObjects() + '\n')
			makefile.write('\t' + libCommand + ' -r ' + v.coinLibrary() + ' '  + v.coinObjects() + '\n')
			if ranlibCommand != '':
				makefile.write('\t' + ranlibCommand + ' ' + v.coinLibrary() + '\n')
			makefile.write('\n')

		makefile.write(v.library() + ': ' + v.objects() + '\n')
		makefile.write('\t' + libCommand + ' -r ' + v.library() + ' '  + v.objects() + ' $(LIBS)\n')
		if ranlibCommand != '':
			makefile.write('\t' + ranlibCommand + ' ' + v.library() + '\n')
		makefile.write('\n')

		makefile.write(v.var + ': ' + v.library())
		if useCoin:
			makefile.write(' ' + v.coinLibrary())
		makefile.write('\n\n')

	makefile.write('clean' + v.var + ':\n')
#	makefile.write('\t' + rmCommand + ' ' + v.objects() + ' ' + v.library() + '\n\n')
	makefile.write('\t' + rmCommand + ' ' + v.call() + '\n\n')

	if useCoin:
		makefile.write('clean' + v.var + '-coin:\n')
		makefile.write('\t' + rmCommand + ' ' + v.coinLibrary() + ' ' + v.coinSharedLibrary() + ' ' + v.call() +'/src/coin\n\n')


	if installPrefix:
		makefile.write('install' + v.var + ':\n')
		makefile.write('\tinstall -d $(DESTDIR)' + installPrefix + '/lib\n')
		if sharedLib:
			makefile.write('\tinstall -m 0644 ' + v.sharedlibrary() + ' $(DESTDIR)' + installPrefix + '/lib/\n')
		else:
			makefile.write('\tinstall -m 0644 ' + v.library() + ' $(DESTDIR)' + installPrefix + '/lib/\n')
		if useCoin:
			if sharedLib:
				makefile.write('\tinstall -m 0644 ' + v.coinSharedLibrary() + ' $(DESTDIR)' + installPrefix + '/lib/\n')
			else:
				makefile.write('\tinstall -m 0644 ' + v.coinLibrary() + ' $(DESTDIR)' + installPrefix + '/lib/\n')
		makefile.write('\n')

def InstallHeaders(curdir, makefile, installPrefix):
	makefile.write('\tinstall -d $(DESTDIR)' + installPrefix + '/' + curdir + '\n')
	names = os.listdir(curdir)
	for name in names:
		filename = curdir + '/' + name
		if (os.path.isdir(filename)):
			InstallHeaders(filename, makefile, installPrefix)
		elif name.endswith('.h') or name.endswith('.hpp'):
			makefile.write('\tinstall -m 0644 ' + filename + ' $(DESTDIR)' + installPrefix + '/' + curdir + '\n')

if installPrefix:
	makefile.write('install: installrelease install-headers install-pkgconfig\n\n')
	makefile.write('install-headers:\n')
	InstallHeaders('include/ogdf', makefile, installPrefix)
	if useCoin:
		InstallHeaders('include/coin', makefile, installPrefix)
	makefile.write('\ninstall-pkgconfig: ogdf.pc\n')
	makefile.write('\tinstall -d $(DESTDIR)' + installPrefix + '/lib/pkgconfig\n')
	makefile.write('\tinstall -m 0644 ogdf.pc $(DESTDIR)' + installPrefix + '/lib/pkgconfig\n')

makefile.write('\ndistclean: clean-doc')
for v in versions:
	makefile.write(' clean' + v.var)
makefile.write('\n\trm -rf Makefile ogdf.pc include/ogdf/internal/config_autogen.h include/coin/config.h\n')

makefile.close()

print('Makefile generated')


