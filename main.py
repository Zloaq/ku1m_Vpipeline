#!/Users/motomo/opt/anaconda3/bin/python3

import os
import sys
import re
import glob
import shutil
import subprocess
import signal
import numpy as np
from astropy.io import fits
from pyraf import iraf

import download
import bottom
import lcs
import dostack_auto as dos_a
import comb_phot as com_p


class readparam():

	def __init__(self, date, objectname):
		
		path_program = os.path.abspath(__file__)
		dir_of_program = os.path.dirname(path_program)
		dir1 = os.path.join(dir_of_program, 'img_proce.param')
		if not os.access(dir1, os.R_OK):
			print('img_proce.param is not found')
			sys.exit()

		with open(dir1) as f:
			for line in f:
				if line.startswith('#'):
					continue
				varr = line.split()
				items = ''
				if varr[1].startswith('\'') or varr[1].startswith('\"'):
#					varr[1] = varr[1].replace('\'', '', 1)
					for index, item in enumerate(varr[1:]):
						if item.endswith('\'') or item.endswith('\"'):
#							item = item.replace('\'', '', 1)
							items = items + item
							varr[1] = items
							break
						items = items + item + ' '
				if '{' in varr[1]:
					varr[1] = varr[1].format(date = date, objectname = objectname)
				exec("self." + varr[0] + " = " + varr[1])
				


class readobjfile():
	
	def __init__(self, param, objectname):
		
		jhkcoordslist = []
		gicoordslist = []

		dir1 = os.path.join(param.objname_dir, objectname)
		if not os.access(dir1, os.R_OK):
			print(objectname, 'file is not found')
			sys.exit()

		with open(dir1, 'r') as f:
			section = None
			for index, line in enumerate(f):
				linelist = line.split()
				if line.startswith('#') or not linelist:
					continue
				if line.startswith('[') and line.endswith(']'):
					section = line[1:-1]
					continue
				if linelist[0] == 'SearchName':
					self.SearchName = []
					for parampart in linelist[1:]:
						if parampart.startswith('#'):
							break
						exec("self." + linelist[0] + ".append(" + parampart + ")")
					continue

				if len(linelist) == 1 or line[1].startswith('#'):
					continue
				elif section == 'Parameters':
					exec("self." + linelist[0] + " = " + linelist[1])
			
				if section == 'Coords_jhk':
					jhkcoordslist.append(f)

				if section == 'Coords_gi':
					gicoordslist.append(f)
		
		if jhkcoordslist:
			dir2 = os.path.join(param.work_dir, 'Reference_jhk.coo')
			with open(dir2, 'w') as f2:
				f2.write('#\n')
				for coo in jhkcoordslist:
					f2.write(coo)
					f2.write('\n')

		if gicoordslist:
			dir3 = os.path.join(param.work_dir, 'Reference_gi.coo')
			with open(dir3, 'w') as f3:
				f3.write('#\n')
				for coo in gicoordslist:
					f3.write(coo)
					f3.write('\n')

'''
				paramlist = []
				for parampart in linelist[1:]:
					if parampart.startswith('#'):
						break
					paramlist.append(parampart)

				if len(paramlist) == 0:
					continue
				if linelist[0] == 'SearchName':
					self.linelist[0] = []
					for parampart in paramlist:
						exec("self." + linelist[0] + ".append(" + parampart + ")")

				elif len(paramlist) == 1:
					exec("self." + linelist[0] + " = " + paramlist[0])
'''


class readlog():
	
	def __init__(self, filename):
		
		if not os.access(filename, os.R_OK):
			print(filename,' is not found')
			self.log = 'nothing'
		
		else:
			self.log = 'exists'
			with open(filename) as f:
				for line in f:
					if line.startswith('#'):
						continue
					varr = line.split()
					if len(varr) <= 1 :
						continue
					exec("self." + varr[0] + " = " + varr[1])



class readheader():
    
    def __init__(self, fitslist):
        
        self.object = []
        self.exptime = []
        self.jd = []
        self.mjd = []
        self.ra = []
        self.dec = []
        self.airmass = []
        self.offsetra = []
        self.offsetde = []
        self.offsetro = []
        self.focus = []
        self.azimuth = []
        self.altitude = []
        self.rotator = []
        
        for f1 in fitslist:
            
            hdu = fits.open(f1)
            self.object.append(hdu[0].header['OBJECT'])
            self.exptime.append(hdu[0].header['EXPTIME'])
            self.jd.append(hdu[0].header['JD'])
            self.mjd.append(hdu[0].header['MJD'])
            self.ra.append(hdu[0].header['RA'])
            self.dec.append(hdu[0].header['DEC'])
            self.airmass.append(hdu[0].header['AIRMASS'])
            self.offsetra.append(hdu[0].header['OFFSETRA'])
            self.offsetde.append(hdu[0].header['OFFSETDE'])
            self.offsetro.append(hdu[0].header['OFFSETRO'])
            self.focus.append(hdu[0].header['FOCUS'])
            self.azimuth.append(hdu[0].header['AZIMUTH'])
            self.altitude.append(hdu[0].header['ALTITUDE'])
            self.rotator.append(hdu[0].header['ROTATOR'])
            hdu.close()


def match_object(fitslist, search_name_list):
	
	obnamelist = []
	wrong_index = []
	for index, fits_file in enumerate(fitslist):
		try:
			HDUlist = fits.open(fits_file)
			obnamelist.append(HDUlist[0].header['OBJECT'])
		except:
			print('header is broken', fits_file)
			wrong_index.append(index)
	
	for index in wrong_index[::-1]:
		del fitslist[index]
		del obnamelist[index]
		
		
	outfitslist = []
	outobnamelist = []
	found_index = []
	for name1 in search_name_list:
		for index, name2 in enumerate(obnamelist):
			if name1 in name2:
				outfitslist.append(fitslist[index])
				outobnamelist.append(obnamelist[index])
				found_index.append(index)
		for index in found_index[::-1]:
			del fitslist[index]
			del obnamelist[index]
		found_index.clear()
			
	return outfitslist, outobnamelist
	

def glob_latestproc(bands, fitspro):
	
	fitslist = []
	key = '.fits'
	for i in fitspro[::-1]:
		key = '_' + i + key
	for band in bands:
		search_pattern = f'{band}*{key}'
		fitslist += glob.glob(search_pattern)
	fitslist.sort()

	return fitslist


def add_latestproc(inlist):
	
	outlist = []
	for i in inlist():
		i.split('_')
		outlist.append(re.sub(r'.fits'))

	

def write_log(filename, paramdict):
    
	# fitspro は中で展開可
	# paramdict = [param_name:param(str or list)]
	dict1 = {}
	linelist2 = []
	if os.access(filename, os.R_OK):
		with open(filename, 'r') as logr:
			linelist = logr.readlines()
		for line in linelist:
			if line.startswith('#') or not line.strip():
				linelist2.append(line)
			else:
				line0 = line.split('#', 1)
				list1 = line0[0].split()
				dict1[list1[0]] = list1[1]
				linelist2.append([list1[0], dict1[list1[0]], line0[1:]])

	for key in paramdict:
		if key in dict1:
			dict1[key] = paramdict[key]
		else:
			linelist2.append([key, paramdict[key]])

	with open(filename, 'w') as logw:
		for line in linelist2:
			if line.startswith('#') or not line.strip():
				logw.write(line + '\n')
			else:
				param_name = '{:<16}'.format(line[0])
				logw.write(param_name)
				if isinstance(line[1]):
					param_part = ''
					for varr in line[1]:
						param_part = param_part + varr
					else:
						param_part = line[1]
				logw.write(param_part)
				if len(line) > 3:
					logw.write('   #')
					logw.write(line[2:])
				logw.write('\n')



def make_objfile(objectname, param):

	objfile = os.path.join(param.objname_dir, objectname)

	if os.access(objfile, os.R_OK):
		print(f'{objectname} file is already exists.')
		sys.exit()

	with open(objectname, 'w') as file:
		file.write('\n')
		file.write('[Parameters]\n')
		file.write('# Please add information of OBJECT.' + '\n')
		file.write('# Please separate names with a space.' + '\n')
		file.write('#  ex) \'Searchname1\' \'Searchname2\' \'Searchname3\'' + '\n\n')
		file.write('{:<16}'.format('ObjectName') + f'\'{objectname}\'' + '\n')
		file.write('{:<16}'.format('SearchName') + f'\'{objectname}\'' +'\n')
		file.write('\n\n')
		file.write('[Coords_jhk]\n')
		file.write('# Please add reference coordinates (in iraf format) for aperture phot.\n')
		file.write('# The order of coordinates is Obj → C1 → C2 → others.\n')
		file.write('# Please list as many coordinates as possible.\n')
		file.write('\n')
		file.write('[Coords_gi]\n')
		file.write('# Please add reference coordinates (in iraf format) for aperture phot.\n')
		file.write('# The order of coordinates is Obj → C1 → C2 → others.\n')
		file.write('# Please list as many coordinates as possible.\n')
		file.write('\n')
	
	print(f'{objectname} file is created')
		


def execute_code(param, objparam, log, bands='jhk'):
	
	bottom.setparam()
	if log.log == 'exits':
		fitspro = log.fitspro
	else:
		fitspro = []
		subprocess.run('rm *.fits', shell=True)
		iraf.chdir(param.rawdata_dir)
		fitslist = glob.glob('*.fits')
		fitslist, obnamelist = match_object(fitslist, objparam.SearchName)
		for file_name in fitslist:
			shutil.copy(file_name, param.work_dir)
	
	iraf.chdir(param.work_dir)

	if param.quicklook == 1:
		print('yet yet yet')
		sys.exit()

	if param.flatdiv == 1:
		print('yetyetyet')
		sys.exit()
		fitslist = glob_latestproc(bands, fitspro)
		lcs.flat_division(fitslist)
		fitspro.append('fl')
	
	if param.cut == 1:
		fitslist = glob_latestproc(bands, fitspro)
		bottom.cut_copy(fitslist, '[9:328,9:264]', './')
		fitspro.append('cut')

	if param.skylevsub == 1:
		fitslist = glob_latestproc(bands, fitspro)
		lcs.method2_1(fitslist)
		fitspro.append('lev')
	elif param.yskylevsub == 1:
		fitslist = glob_latestproc(bands, fitspro)
		lcs.method2_2(fitslist)
		fitspro.append('ylev')
	elif param.skylevsub == 1 and param.yskylevsub == 1:
		fitslist = glob_latestproc(bands, fitspro)
		lcs.method2_1(fitslist)
		fitspro.append('lev')
		
	if param.skysub == 1:
		fitslist = glob_latestproc(bands, fitspro)
		header = readheader(fitslist)
		fitslist = lcs.method3(fitslist, header.object)
		header = readheader(fitslist)
		lcs.method4(fitslist, header.object)
		fitspro.append('sky')
		
	if param.dostack == 1:
		fitslist = glob_latestproc(param.dostack_base, fitspro)
		header = readheader(fitslist)
		dos_a.dostack_main(fitslist, param)
		fitspro.append('geo')


	if param.comb_per_set == 1:
		print('combs')
		fitslist = glob_latestproc(bands, fitspro)
		com_p.comb_pset(fitslist, argvs[1], argvs[2])

	if param.comb_all == 1:
		print('comb_all')
		fitslist = glob_latestproc(bands, fitspro)
		com_p.comb_all(fitslist, argvs[1], argvs[2])

	if param.aperture_phot == 1:
		print('phot')
		


if __name__ == '__main__':
	
	argvs = sys.argv
	argc = len(argvs)
	fitspro = []
	
	print(argc)

	if argc == 2:
		#object ファイルの生成
		param = readparam('', argvs[1])
		os.makedirs(param.objname_dir, exist_ok=True)
		iraf.chdir(param.objname_dir)
		make_objfile(argvs[1], param)
	
	elif argc == 3:
		param = readparam(argvs[1], argvs[2])
		objparam = readobjfile(param, argvs[2])

		path = os.path.join(param.work_dir)
		os.makedirs(path, exist_ok=True)
		iraf.chdir(path)
		log = readlog('log.txt')
		execute_code(param, objparam, log)

	elif argc == 4:
		print('yetyetyet, only jhk')
		sys.exit()
		param = readparam(argvs[1], argvs[2])
		objparam = readobjfile(param, argvs[2])

		path = os.path.join(param.work_dir)
		os.makedirs(path, exist_ok=True)
		iraf.chdir(path)
		log = readlog('log.txt')
		execute_code(param, objparam, log, argvs[3])
	
	else:
		print('usage1 ./main.py [OBJECT file name] ')
		print('usage2 ./main.py [YYMMDD][object name]')