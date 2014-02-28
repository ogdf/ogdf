#!/usr/bin/env python
# Make Makefile
#
# October 2012
# Markus Chimani, markus.chimani@cs.tu-dortmund.de
# Carsten Gutwenger, carsten.gutwenger@cs.tu-dortmund.de
#########################################################

# use Python 3 print function on both Python 2.x and Python 3.x
from __future__ import print_function

import os, sys, fnmatch, posixpath, shutil, subprocess

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

def loadConfig(sect, key, noError = False):
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
	config.readfp(open('makeMakefile.config'))
except IOError:
	bailout('makeMakefile.config not found')

if not config.has_section('GENERAL'):
	bailout('Section "GENERAL" is missing')
if not config.has_section('OGDF'):
	bailout('Section "OGDF" is missing')
if not config.has_section('VERSIONS'):
	bailout('Section "VERSIONS" is missing')
if not config.has_section('COIN'):
	bailout('Section "COIN" is missing')

compilerCommand = loadConfig('GENERAL', 'compilerCommand')
compilerFlags = loadConfig('GENERAL', 'compilerParams')
libCommand = loadConfig('GENERAL', 'libCommand')
sharedlibCommand = loadConfig('GENERAL', 'sharedlibCommand')
rmCommand = loadConfig('GENERAL', 'rmCommand')
mkdirCommand = loadConfig('GENERAL', 'mkdirCommand')
ranlibCommand = loadConfig('GENERAL', 'ranlibCommand', True)
installPrefix = loadConfig('GENERAL', 'installPrefix', True)

sharedLib = loadConfig('GENERAL', 'sharedLib').startswith('t')
libName = loadConfig('OGDF', 'libName')
sharedlibName = loadConfig('OGDF', 'sharedlibName')
includeLegacyCode = loadConfig('OGDF', 'includeLegacyCode').startswith('t')
memoryManager = loadConfig('OGDF', 'memoryManager', True)

# defines for auto generated config_autogen.h header file
config_defines = ''

gccMessageLength = loadConfig('GENERAL', 'gccMessageLength', True)
if gccMessageLength:
	compilerFlags = '-fmessage-length=' + gccMessageLength + ' ' + compilerFlags

libs = ''

if sharedLib:
	config_defines += '#define OGDF_DLL\n'
	if sys.platform == 'win32' or sys.platform == 'cygwin':
		libs += ' -lpsapi'
	else:
		compilerFlags += ' -fPIC'

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
	solverLDFlags   = loadConfig('COIN', 'solverLDFlags')

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
		libs += ' -lCOIN'

ogdfFlags = '-I./include ' + addIncludes
coinFlags = '$(COIN_INSTALL_DEFINES) -I./include/coin ' + addIncludes

if sharedLib:
	ogdfFlags += ' -DOGDF_INSTALL'

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

print('Resulting compiler call (OGDF):', ' '.join([compilerCommand, compilerFlags, ogdfFlags]))
if useCoin:
	print('Resulting compiler call (COIN):', ' '.join([compilerCommand, compilerFlags, coinFlags]))

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

# define general variables

makefile.write('CC = ' + compilerCommand + '\n')
makefile.write('CXXFLAGS = ' + compilerFlags + '\n')
makefile.write('OGDFFLAGS = ' + ogdfFlags + '\n')
makefile.write('COINFLAGS = ' + coinFlags + '\n')
if ranlibCommand:
	makefile.write('RANLIB = ' + ranlibCommand + '\n')
makefile.write('RM = ' + rmCommand + '\n')
makefile.write('AR = ' + libCommand + '\n')
makefile.write('LD = ' + sharedlibCommand + '\n')
makefile.write('MKDIR = ' + mkdirCommand + '\n')

# define release & debug

for v in versions:
	makefile.write(v.var + ' = ' + v.path() + '\n')
	makefile.write('CXXFLAGS_' + v.var + ' = ' + v.params + '\n')
makefile.write('\n');

# just the def. nothing happens yet
def Walk(curdir):
	objs = []
	names = os.listdir(curdir)
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
				objs = objs + Walk(fullname)
		else:
			for pat in [ '*.c', '*.cpp' ]:
				if fnmatch.fnmatch(name, pat):
					objfullname = fullname[:-len(pat)+2] + 'o'
					objs.append(objfullname)

					targetAndDepend = subprocess.check_output(' '.join([compilerCommand, compilerFlags, ogdfFlags, '-MM', fullname]), shell=True)
					for v in versions:
						# print target&depend: add full path spec, incl. version & ignore extra line
						path = v.call() + '/' +fullname[:-len(name)]
						makefile.write(path + targetAndDepend[:-1] + '\n')

						# ensure folder
						makefile.write('\t$(MKDIR) ' + v.call() + '/' + fullname[:-len(name)-1] + '\n')
						# what to do: call the compiler
						makefile.write('\t$(CC) $(CXXFLAGS) $(OGDFFLAGS) $(CXXFLAGS_' + v.var + ') -o $@ -c ' + fullname + '\n\n')

					# pattern found: don't try other suffix
					break
	return objs

def WalkCoin(curdir):
	objs = []
	names = os.listdir(curdir)
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
				objs = objs + WalkCoin(fullname)
		else:
			for pat in [ '*.c', '*.cpp' ]:
				if fnmatch.fnmatch(name, pat):
					objfullname = fullname[:-len(pat)+2] + 'o'
					if not name.startswith('unitTest'): #conflict if in library!
						objs.append(objfullname)

					# skip .h dependency check here

					for v in versions:
						path = v.call() + '/' +objfullname
						makefile.write(path + ': ' + fullname + '\n')

						# ensure folder
						makefile.write('\t$(MKDIR) ' + v.call() + '/' + fullname[:-len(name)-1] + '\n')
						# what to do: call the compiler
						makefile.write('\t$(CC) $(CXXFLAGS) $(COINFLAGS) $(CXXFLAGS_' + v.var + ') -o $@ -c ' + fullname + '\n\n')

					# pattern found: don't try other suffix
					break
	return objs

# OGDF
objs = Walk('./src/ogdf')

# List all Objs for use in lib-generation and clear
for v in versions:
	makefile.write(v.objects()[2:-1] + ' = \\\n')
	for o in objs:
		makefile.write('\t' + v.call() + '/' + o + ' \\\n')
	makefile.write('\n')

# Coin
if useCoin:
	objsCoin = WalkCoin('./src/coin')

	# List all Objs for use in lib-generation and clear
	for v in versions:
		makefile.write(v.coinObjects()[2:-1] + ' = \\\n')
		for o in objsCoin:
			makefile.write('\t' + v.call() + '/' + o + ' \\\n')
		makefile.write('\n')

# generate test binary targets

makefile.write('\n#######################################################')
makefile.write('\n# test binary targets\n\n')

objsTest = Walk('./test')

for v in versions:
	makefile.write('test/test-' + v.var + ': ' + ' '.join([v.call() + '/' + x for x in objsTest]) + '\n')
	makefile.write('\t$(CC) -o $@ $^ -pthread -L' + v.call() + ' -lOGDF')
	if useCoin:
		makefile.write(' -lCOIN ' + solverLDFlags)
	makefile.write('\n\n')

# generate alls and cleans etc...

for v in versions:
	makefile.write('\n#######################################################')
	makefile.write('\n# all, clean, etc. for ' + v.var + '\n\n')

	makefile.write('.PHONY: ' + v.var + ' install' + v.var + ' clean' + v.var + ' clean' + v.var + '-coin\n\n')

	if sharedLib:
		if useCoin:
			makefile.write(v.coinSharedLibrary() + ': ' + v.coinObjects() + '\n')
			makefile.write('\t$(LD) -shared -o ' + v.coinSharedLibrary() + ' ' + v.coinObjects() + '\n\n')

		makefile.write(v.sharedlibrary() + ': ' + v.objects() + '\n')
		makefile.write('\t$(LD) -shared -o ' + v.sharedlibrary() + ' ' + v.objects() + ' ' + libs + ' $(LIBS)\n\n')

		makefile.write(v.var + ': ' + v.sharedlibrary())
		if useCoin:
			makefile.write(' ' + v.coinSharedLibrary())
		makefile.write('\n\n')

	else:
		if useCoin:
			makefile.write(v.coinLibrary() + ': ' + v.coinObjects() + '\n')
			makefile.write('\t$(AR) -r ' + v.coinLibrary() + ' '  + v.coinObjects() + '\n')
			if ranlibCommand:
				makefile.write('\t$(RANLIB) ' + v.coinLibrary() + '\n')
			makefile.write('\n')

		makefile.write(v.library() + ': ' + v.objects() + '\n')
		makefile.write('\t$(AR) -r ' + v.library() + ' '  + v.objects() + ' $(LIBS)\n')
		if ranlibCommand:
			makefile.write('\t$(RANLIB) ' + v.library() + '\n')
		makefile.write('\n')

		makefile.write(v.var + ': ' + v.library())
		if useCoin:
			makefile.write(' ' + v.coinLibrary())
		makefile.write(' test/test-' + v.var + '\n\n')

	makefile.write('clean' + v.var + ':\n')
#	makefile.write('\t$(RM) ' + v.objects() + ' ' + v.library() + '\n\n')
	makefile.write('\t$(RM) ' + v.call() + ' test/test-' + v.var + '\n\n')

	if useCoin:
		makefile.write('clean' + v.var + '-coin:\n')
		makefile.write('\t$(RM) ' + v.coinLibrary() + ' ' + v.coinSharedLibrary() + ' ' + v.call() +'/src/coin\n\n')

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
	makefile.write('.PHONY: install-headers install-pkgconfig\n')
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
makefile.write('\n\t$(RM) Makefile ogdf.pc include/ogdf/internal/config_autogen.h include/coin/config.h\n')

makefile.close()

print('Makefile generated')


