#!/Users/motomo/opt/anaconda3/bin/python3

import os
import sys
import re
import glob
import shutil
import subprocess
import statistics
import multiprocessing
import numpy as np
from astropy.io import fits
from time import sleep
from pyraf import iraf

import bottom


class readparam():

	def __init__(self):
		
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
				exec("self." + varr[0] + " = " + varr[1])
		
				
				
class readlog():
	
	def __init__(self, logfile):
		
		if not os.access(logfile, os.R_OK):
			print(logfile,' is not found')
			sys.exit()
		
		with open(logfile) as f:
			for line in f:
				varr = line.split()
				exec("self." + varr[0] + " = " + varr[1])


class readobname():
    
    def __init__(self, inlist):
        
        self.object = []
        
        for f1 in inlist:
            
            hdu = fits.open(f1)
            self.object.append(hdu[0].header['OBJECT'])
            hdu.close()



class readmjd():
    
	def __init__(self, inlist):
		
		self.mjd = []
		for f1 in inlist:
			print(f1)
			hdu = fits.open(f1)
			self.mjd.append(hdu[0].header['MJD'])
			hdu.close()



def opends9():
    
    try:
        subprocess.run('xpaset -p ds9 exit', shell=True)
    except:
        print('ok')
    subprocess.run('ds9&', shell=True)
    sleep(5)



def combine(inlist, skyname, method, zero):
	
    iraf.imcombine.unlearn()
    iraf.imcombine.combine = method
    iraf.imcombine.reject = 'sigclip'
    iraf.imcombine.zero = zero
    iraf.imcombine.mclip = 'yes'
    
    iraf.imcombine('@' + inlist, skyname, Stdout=1)



def imstat(infile, fields):

	iraf.imstatistics.unlearn()
	iraf.imstatistics.fields = fields
	iraf.imstatistics.lower = -1000
	iraf.imstatistics.nclip = 5
	iraf.imstatistics.format = 'No'
	return iraf.imstatistics(infile, Stdout=1)



def daofind(inputf, outputf, sigmav):
    
	iraf.noao()
	iraf.digiphot(Stdout=1)
	iraf.apphot(Stdout=1)
	iraf.apphot.datapars.scale = 0.692
	iraf.apphot.datapars.fwhmpsf = fwhm_med
	iraf.apphot.datapars.sigma = sigmav
	iraf.apphot.datapars.datamin = 0
	iraf.apphot.datapars.datamax = 20000
	iraf.apphot.findpars.threshold = param.threshold
	iraf.daofind.unlearn()
	iraf.daofind.output = outputf
	iraf.daofind.interac = 'no'
	iraf.daofind.verify = 'no'
	iraf.daofind.mode = 'h'
	iraf.daofind(inputf)



def phot(inimage, outfile, coords, aperture, annulus):
	
	iraf.apphot.fitskypars.annulus = annulus
	iraf.apphot.photpars.aperture = aperture
	iraf.apphot.phot.image = inimage
	iraf.apphot.phot.coords = coords
	iraf.apphot.phot.output = outfile
	iraf.apphot.phot.skyfile = ''
	iraf.apphot.phot.interactive = 'No'
	iraf.apphot.phot.verify = 'No'
	iraf.apphot.phot.verbose = 'No'
	iraf.apphot.phot.mode = 'h'
	iraf.apphot.phot(inimage)

	

def setparam():
	
	iraf.noao()
	iraf.digi(Stdout=1)
	iraf.app(Stdout=1)
	iraf.images()
	iraf.immatch()
	iraf.tv()
	iraf.digiphot.unlearn()
	iraf.apphot.unlearn()
	iraf.ptools.unlearn()
	iraf.txdump.unlearn()
	iraf.apphot.datapars.scale = 0.692
	iraf.apphot.datapars.fwhmpsf = fwhm_med
	iraf.apphot.datapars.datamax = 20000
	iraf.apphot.datapars.datamin = 0
	iraf.apphot.datapars.readnoise = 50
	iraf.apphot.datapars.epadu = 4
	
	iraf.apphot.centerpars.calgorithm = 'centroid'
	
	iraf.apphot.fitskypars.dannulus = 2.
	iraf.apphot.findpars.threshold = param.threshold
	iraf.apphot.findpars.sharphi = 0.8
	iraf.apphot.photpars.zmag = 0


def txdump(infile, fields):
	
	iraf.txdump.fields = fields
	iraf.txdump.mode = 'h'
	return iraf.txdump(infile, Stdout=1)
	

def phmethod0_1():
	
	subprocess.run(f'rm {argvs[2]}.coo', shell=True)
	lev = imstat('j'+argvs[1]+'-0000'+ist0.fitspro+'_all.fits', 'stddev')
	print(lev)
	daofind('j'+argvs[1]+'-0000'+ist0.fitspro+'_all.fits', argvs[2]+'.coo', lev[0])


def phmethod0_2():
	
	opends9()
	subprocess.run(f'rm {argvs[2]}.coo', shell=True)
	list1 = iraf.imexam('j'+argvs[1]+'-0000'+ist0.fitspro+'_all.fits', 1, Stdout=1)
	fwhm1 = []
	with open(f'{argvs[2]}.coo' ,'w') as f1:
		for line in list1:
			if line.startswith('#'):
				continue
			prof = line.split()
			if len(prof) == 4:
				f1.write(prof[0]+'\t'+prof[1]+'\n')
			elif len(prof) == 11:
				fwhm1.append(prof[-2])
	
	global fwhm_med		
	fwhm_med = statistics.median(fwhm1)
	setparam()
	
	

def phmethod1():
	
	subprocess.run('rm *mag.1', shell=True)
	list1 = glob.glob('*' + ist0.fitspro + '.fits')
	list1.sort()
	mjd1 = readmjd(list1)
	for i1 in list1:
		print('phot', i1, i1+'.mag.1', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))
		phot(i1, i1+'.mag.1', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))


def phmethod2():
    # matching objectname
	subprocess.run('rm *n*_all.fits', shell=True)
	subprocess.run('rm *mag.2', shell=True)
	ist1 = readobname(fitslist1)
	findex=0
	for i1 in range(len(ist1.object)):
		if i1+1 == len(ist1.object) or ist1.object[i1] != ist1.object[i1+1]:
			with open('_inlist', 'w') as f:
				for i2 in range(findex, i1+1):
					f.write(fitslist1[i2] + '\n')
			band = fitslist1[findex][1]
			out1 = band + '_' + ist1.object[findex] + '_all.fits'
			combine('_inlist', out1, 'average', 'none')
			findex = i1+1
	list1 = glob.glob('*n*_all.fits')
	for i1 in list1:
		list1.sort()
		print('phot', i1, i1+'.mag.2', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))
		try:
			phot(i1, i1+'.mag.2', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))
		except:
			print('skip')
		 
def phmethod3():
	
	subprocess.run('rm *mag.3', shell=True)
	list1 = glob.glob('*-0000*_all.fits')
	for i1 in list1:
		list1.sort()
		print('phot', i1, i1 + '.mag.3', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))
		phot(i1, i1 + '.mag.3', argvs[2]+'.coo', 2.5*float(fwhm_med), 2.5*float(fwhm_med))



def phot_main():
	
	if param.imexcoords == 1:
		phmethod0_2()
	
	elif param.daocoords == 1:
		phmethod0_1()
	
	if param.eachphot == 1:
		phmethod1()
		
	if param.setcombphot == 1:
		phmethod2()
		
	if param.allcombphot == 1:
		phmethod3()
		
	
def light_curve():
	
	if param.eachphot == 1:
		
		for band in 'jhk':
			list1 = glob.glob(band+'*.mag.1')
			list1.sort()
			inst = readmjd([element[:-6] for element in list1])
			with open(band+'eachphot', 'w') as f1:
				list3 = {}
				for i1 in range(len(list1)):
					print(list1[i1])
					list2 = iraf.txdump(list1[i1], 'id, mag, merr', Stdout=1)
					for i2 in list2:
						i3 = i2.split()
						list3[i3[0]] = i3[1]+' '+i3[2]+' '
					for i4 in [1, 2, 3]:
						f1.write(list3[i4])
					f1.write(f'{inst.mjd[i1]}\n')
					
	if param.setcombphot == 1:
		
		for band in 'jhk':
			list1 = glob.glob(band+'*.mag.2')
			list1.sort()
			inst = readmjd(list1[0:][:-6])
			with open(band+'setcombphot', 'w') as f1:
				list3 = {}
				for i1 in range(len(list1)):
					list2 = txdump(list1[i1], 'id, mag, merr', Stdout=1)
					for i2 in list2:
						i3 = i2.split()
						list3[i3[0]] = i3[1]+' '+i3[2]+' '
					for i4 in [1, 2, 3]:
						f1.write(list3[i4])
					f1.write(f'{inst.mjd[i1]}\n')

		
	if param.allcombphot == 1:
		
		for band in 'jhk':
			list1 = glob.glob(band+'*.mag.3')
			list1.sort()
			inst = readmjd(list1[0:][:-6])
			with open(band+'allcombphot', 'w') as f1:
				list3 = {}
				for i1 in range(len(list1)):
					list2 = txdump(list1[i1], 'id, mag, merr', Stdout=1)
					for i2 in list2:
						i3 = i2.split()
						list3[i3[0]] = i3[1]+' '+i3[2]+' '
					for i4 in [1, 2, 3]:
						f1.write(list3[i4])
					f1.write(f'{inst.mjd[i1]}\n')
	
	
	gnuscript = f"""
		set terminal png
		set output 'eachphot.png'
		set autoscale 
		set grid
		plot "jeachphot" using 7:($1-$3):(sqrt($2*2+$4*2)) with yeerrorbars lc 3\
			 "heachphot" using 7:($1-$3):(sqrt($2*2+$4*2)) with yeerrorbars lc 2\
			 "keachphot" using 7:($1-$3):(sqrt($2*2+$4*2)) with yeerrorbars lc 1
		"""
	print('gnuplot run')
	subprocess.run(["gnuplot"], input=gnuscript, text=True)

	gnuscript2 = f"""
		set terminal png
		set output 'eachphot_2.png'
		set autoscale 
		set grid
		plot "jeachphot" using 7:($5-$3):(sqrt($6*2+$4*2)) with yeerrorbars lc 3\
			 "heachphot" using 7:($5-$3):(sqrt($6*2+$4*2)) with yeerrorbars lc 2\
			 "keachphot" using 7:($5-$3):(sqrt($6*2+$4*2)) with yeerrorbars lc 1
		"""
	print('gnuplot run')
	subprocess.run(["gnuplot"], input=gnuscript2, text=True)



if __name__ == "__main__":

	argvs = sys.argv
	argc = len(argvs)
	param = readparam()
	pattern1 = os.path.join(param.date_dir, argvs[1])
	pattern2 = os.path.join(param.date_dir, argvs[1], argvs[2])
	iraf.chdir(pattern2)
	ist0 = readlog(argvs[2]+'.param')
	fitslist1 = glob.glob('*'+ist0.fitspro+'.fits')
	fitslist1.sort()

	if argc == 4:
		phot_main()
	elif argc == 3:
		light_curve()
	else:
		print('aperture_phot.py [DAY] [obname] [band]')
		print('usage: aperture_phot.py 230103 ZTF23aaaatwl jhk')
		sys.exit()