#!/Users/motomo/opt/anaconda3/bin/python3

import os
import sys
import re
import glob
import fnmatch
import shutil
import signal
import subprocess
import statistics
import multiprocessing
import time
import numpy as np
from astropy.io import fits
from pyraf import iraf

import bottom


class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
	iraf.flpr()
	bottom.setparam()
	raise TimeoutError
signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(0)


def daofind_list(inlist, outlist, fwhms, threshold):
	
	if isinstance(fwhms, list):

		if len(fwhms) != len(inlist):
			print('len(fwhms) != len(inlist')
			sys.exit()

		for i in range(len(inlist)):
			if fwhms[i] == 0:
				continue
			print(len(inlist), len(outlist))
			print(inlist[i], outlist[i], fwhms, threshold)
			bottom.daofind(inlist[i], outlist[i], fwhms[i], threshold)
		
	else:
		for i in range(len(inlist)):
			bottom.daofind(inlist[i], outlist[i], fwhms, threshold)



def daofind_result(daolist, fitslist, nxblock, nyblock):
	
	header = bottom.readheader(fitslist)
	
	dictlist = []
	for index in range(len(daolist)):
		with open(daolist[index], 'r') as f1:
			lines = f1.readlines()
		result = {}
		result['fitsname'] = fitslist[index]
		result['xcen'] = []
		result['ycen'] = []
		result['mag']  = []

		if float(header.offsetra[index]) > 0:
			xmin = float(header.offsetra[index])/0.692 + 5
			xmax = nxblock
		else:
			xmin = 0
			xmax = int(nxblock) - float(header.offsetra[index])/0.692 -5

		if float(header.offsetde[index]) < 0:
			ymin = -float(header.offsetra[index])/0.692 + 5
			ymax = nyblock
		else:
			ymin = 0
			ymax = int(nyblock) - float(header.offsetde[index])/0.692 -5

		for line in lines:
			line_list = line.strip().split()
			if line.startswith('#'):
				continue

			xcen, ycen = float(line_list[0]), float(line_list[1])
			
			if len(line_list) == 8 and xmin < xcen < xmax  and  ymin < ycen < ymax:
				result['xcen'].append(xcen)
				result['ycen'].append(ycen)
				result['mag'].append(float(line_list[2]))
		
		if result['mag']:
			# fwhm measure で25パーセンタイルをとるために並べ替える。
			indices = sorted(range(len(result['mag'])), key=lambda i: result['mag'][i], reverse=False)
			result['xcen'] = [result['xcen'][i] for i in indices]
			result['ycen'] = [result['ycen'][i] for i in indices]
			result['mag']  = [result['mag'][i] for i in indices]
            

		dictlist.append(result)
		
	return dictlist



def fwhmmeasure(dictlist, nxblock, nyblock):
	
	original_stdgraph = iraf.envget("stdgraph") #2024/05/14
	iraf.set(stdgraph="eps") #2024/05/14
	epslist1 = glob.glob('*.eps')

	for dict1 in dictlist:
		if not dict1['xcen']:
			dict1['fwhm'] = 3
			print('block',dict1)
			continue
#		print('pass',dict1)
		fitsname = dict1['fitsname']
		band = fitsname[0]
		psfs = []
		fwhms = []
		
		for index in range(len(dict1['xcen'])):
			bottom.imshift(
			fitsname, 'temp.fits',
			float(nxblock)/2-dict1['xcen'][index], float(nyblock)/2-dict1['ycen'][index]
			)
#			print('psfme start')
			#iraf.display('temp.fits', 1)
			#time.sleep(2)
#			print('\n')
			print(dict1['fitsname'], dict1['xcen'][index], dict1['ycen'][index])
#			print('\n')
			try:
				signal.alarm(1)
				psfresult = bottom.psfmeasure('temp.fits')
#				print('result    ', psfresult)
				print('measure ', fitsname)
				signal.alarm(0)
			except TimeoutError:
				print('Timeouttttttttttttttttt')
				signal.alarm(0)
				os.remove('temp.fits')
				bottom.setparam(band)
				continue
			signal.alarm(0)
			psfs.append(psfresult)
			print('psfme end', dict1['fitsname'])
			os.remove('temp.fits')
		for psf in psfs:
			for line in psf:
				line_list = line.split()
				if not line_list:
					continue
				if line_list[0] == 'Full':
					fwhms.append(float(line_list[-1]))
		
		if fwhms:
			size = len(fwhms)
			sorted_fwhms = sorted(fwhms)
			index = int(size * 0.25)
			dict1['fwhm'] = fwhms[index]
	

	iraf.set(stdgraph=original_stdgraph) #2024/05/14
	epslist2 = glob.glob('*.eps')
	new_eps = list(set(epslist2) - set(epslist1))
#	for file_name in new_eps:
	for file_name in epslist2:
		if os.path.exists(file_name):
			os.remove(file_name)
#			print(f"{file_name} を削除しました。")

	return dictlist



def geotparam(file_list):

	param_list = []

	for filename in file_list:
		geotp = {}
		geotp['fitsid']=filename[1:-4]

		with open(filename, 'r') as f1:

			for line in f1:
				
				line_list = line.split()

				if 'xrefmean' in line:
					geotp['xrefmean']=float(line_list[1])
				elif 'yrefmean' in line:
					geotp['yrefmean']=float(line_list[1])
				elif 'xmean' in line:
					geotp['xmean']=float(line_list[1])
				elif 'ymean' in line:
					geotp['ymean']=float(line_list[1])
				elif 'xrotation' in line:
					geotp['xrotation']=float(line_list[1])
				elif 'yrotation' in line:
					geotp['yrotation']=float(line_list[1])
				elif 'xshift' in line:
					geotp['xshift']=float(line_list[1])
				elif 'yshift' in line:
					geotp['yshift']=float(line_list[1])
		
		param_list.append(geotp)
					
	return param_list


def domethod1(inlist, fwhms, threshold, param):
	#半値幅を検出の数が多いものでざっくり
	outlist = {}
	results = {}
	for fwhm in fwhms:
		outlist[fwhm] = []
		for i2 in inlist:
			outlist[fwhm].append(i2[:14]+'_' +str(fwhm)+'.coo')
			
	stddevlist1 = bottom.imstat(inlist, 'stddev')
	hdulist = [fits.getheader(file, ext=0) for file in inlist]
	medlist = [float(hdu['MEDIAN']) for hdu in hdulist]
	
	if param.multi == 1:
		process1 = []
		for fwhm in fwhms:
			p = multiprocessing.Process(
			target=daofind_list, args=(inlist, outlist[fwhm], fwhm, threshold)
			)
			process1.append(p)
			p.start()
		for p in process1:
			p.join()

	elif param.multi == 0:
		for fwhm in fwhms:
			daofind_list(inlist, outlist[fwhm], fwhm, threshold)
	
	xcen_counter = {}
	
	data =  fits.getdata(inlist[0])
	nxblock = len(data[0])
	nyblock = len(data)
	
	for fwhm in fwhms:
		
		results[fwhm] = daofind_result(outlist[fwhm], inlist, nxblock, nyblock)
		xcen_counter[fwhm] = [len(result['xcen']) for result in results[fwhm]]
		
#	goodfwhm = max(xcen_counter, key=lambda k: statistics.median(xcen_counter[k]))
	goodfwhm = min(xcen_counter, key=lambda k: statistics.stdev(xcen_counter[k]))

	return results[goodfwhm], goodfwhm


def domethod2(indictlist, firstfwhm, roop, threshold):
	# fwhm イテレーション
	fwhm = firstfwhm
	for num in range(roop):
		data =  fits.getdata(indictlist[0]['fitsname'])
		nxblock = len(data[0])
		nyblock = len(data)
		fitslist = [dict['fitsname'] for dict in indictlist]
		outlist = [re.sub(r'.fits', r'.coo', dict['fitsname']) for dict in indictlist]
		daofind_list(fitslist, outlist, fwhm, threshold)
		indictlist = daofind_result(outlist, fitslist,  nxblock, nyblock)
		indictlist = fwhmmeasure(indictlist, nxblock, nyblock)
		fwhm = [dict['fwhm'] for dict in indictlist]

	return indictlist, fwhm


def domethod3_1(coordslist, bands):

	for filename in coordslist[1:]:
		print('xyxymatch', filename, coordslist[0], re.sub(r'.coo', r'.match', filename))
		iraf.xyxymatch(filename, coordslist[0], re.sub(r'.coo', r'.match', filename), 3.0, Stdout=1)

	listmatch = glob.glob('*.match')
	listmatch.sort()

	data = fits.getdata(re.sub(r'.coo', r'.fits', coordslist[0]))
	nxblock = len(data[0])
	nyblock = len(data)
	for filename in listmatch:
		print('geomap ', filename)
		bottom.geomap(filename, nxblock, nyblock, re.sub(r'.match', r'.geo', filename), 'rotate')
		
	listgeo = glob.glob('*.geo')
	listgeo.sort()

	geotlist = geotparam(listgeo)

	bottom.geotran_param(nxblock, nyblock)
	print('len(geotlist) == ', len(geotlist))
	
	xmean0 = {}
	ymean0 = {}
	xrot = {}
	yrot = {}

	xmean0['j'], ymean0['j'], xrot['j'], yrot['j'] = 0.0, 0.0, 0.0, 0.0
	xmean0['h'], ymean0['h'], xrot['h'], yrot['h'] = -1.55189, -14.77146, 0.13356, 0.13356
	xmean0['k'], ymean0['k'], xrot['k'], yrot['k'] = -19.97968, -9.8439, 359.96445, 359.96445

	xmean0['g'], ymean0['g'], xrot['g'], yrot['g'] = -1.55189, -14.77146, 0.13356, 0.13356
	xmean0['i'], ymean0['i'], xrot['i'], yrot['i'] = -19.97968, -9.8439, 359.96445, 359.96445

#	for dict in geotlist:
	for band in bands:
#		if len(dict) == 1:
#			continue

#		if  10 < dict['xrotation'] < 350 or 10 < dict['yrotation'] < 350:
#			continue

#		for band in bandlist:
		for dict in geotlist:
			if len(dict) == 1:
				print('reject', band + dict['fitsid']+'.fits')
				continue
			if  10 < dict['xrotation'] < 350 or 10 < dict['yrotation'] < 350:
				continue
			iraf.geotran.xout = dict['xrefmean']
			iraf.geotran.yout = dict['yrefmean']
			iraf.geotran.xin = dict['xmean'] + xmean0[band]
			iraf.geotran.yin = dict['ymean'] + ymean0[band]
			iraf.geotran.xrotation = dict['xrotation'] + xrot[band]
			iraf.geotran.yrotation = dict['yrotation'] + yrot[band]
			print('geotran ',band + dict['fitsid'] + '_geo.fits')
			try:
				iraf.geotran(band + dict['fitsid']+'.fits',band + dict['fitsid']+'_geo.fits')
			except:
				print(f"Cannot open image ({band}{dict['fitsid']}.fits)")
	


def domethod3_2(coordslist, bands):

	iraf.geomap.maxiter = 10
	iraf.geomap.reject = 1.5

	for filename in coordslist[1:]:
		print('lnum_match', filename, coordslist[0], re.sub(r'.coo', r'.match', filename))
		bottom.lnum_match(filename, coordslist[0], re.sub(r'.coo', r'.match', filename))

	listmatch = glob.glob('*.match')
	listmatch.sort()

	data = fits.getdata(re.sub(r'.coo', r'.fits', coordslist[0]))
	nxblock = len(data[0])
	nyblock = len(data)
	for filename in listmatch:
		print('geomap ', filename)
		bottom.geomap(filename, nxblock, nyblock, re.sub(r'.match', r'.geo', filename), 'shift')
		
	listgeo = glob.glob('*.geo')
	listgeo.sort()

	geotlist = geotparam(listgeo)

	bottom.geotran_param(nxblock, nyblock)
	print('len(geotlist) == ', len(geotlist))
	
	xmean0 = {}
	ymean0 = {}
	xrot = {}
	yrot = {}

	xmean0['j'], ymean0['j'] = 0.0, 0.0
	xmean0['h'], ymean0['h'] = -1.55189, -14.77146
	xmean0['k'], ymean0['k'] = -19.97968, -9.8439

	for band in bands:
#	for dict in geotlist:
#		for band in bandlist:
		for dict in geotlist:
			try:
				iraf.geotran.xshift = -1 * dict['xshift'] + xmean0[band]
				iraf.geotran.yshift = -1 * dict['yshift'] + ymean0[band]
				iraf.geotran(band + dict['fitsid']+'.fits',band + dict['fitsid']+'_geo.fits')
				print('geotran ',band + dict['fitsid'] + '_geo.fits')
			except:
				print(f"Cannot do geotran of image ({band}{dict['fitsid']}.fits)")
	
	

def dostack_main(fitslist, param):

	subprocess.run('rm *.match' ,shell=True)
	subprocess.run('rm *geo*' ,shell=True)
	subprocess.run('rm *.coo' ,shell=True)
	subprocess.run('rm *_xysh*' ,shell=True)

	dictlist1, fwhm = domethod1(fitslist, [2, 3, 4], 5, param)

	subprocess.run('rm *.coo' ,shell=True)

	dictlist2, fwhmlist = domethod2(dictlist1, fwhm, 1, 5)

	coordscounter = statistics.mean([len(dict['xcen']) for dict in dictlist2])

	if coordscounter < 1:
		print('no star')
		sys.exit()

	if coordscounter < 5:
		print('num of star < 3')
		param.georotate == 0
		param.xyshift == 1


	coordslist = glob.glob('*.coo')
	coordslist.sort()

	if param.georotate == 1:
		domethod3_1(coordslist, 'jhk')

	if param.xyshift == 1:
		domethod3_2(coordslist, 'jhk')
	



if __name__ == "__main__":

	print('This is module.')