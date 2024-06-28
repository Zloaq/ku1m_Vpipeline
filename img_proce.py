#!/usr/bin/env python3

import os
import sys
import re
import glob
import shutil
import warnings
import numpy as np
from astropy.io import fits
from pyraf import iraf


class read_param:
	
	def __init__(self):

		path_program = os.path.abspath(__file__)
		dir_of_program = os.path.dirname(path_program)
		dir1 = os.path.join(dir_of_program, 'img_proce.param')
		if not os.access(dir1, os.R_OK):
			print('img_proce.param is not found')
			sys.exit()
	
		with open(dir1) as f:
			for line in f:
				varr = line.split()
				exec("self." + varr[0] + " = " + varr[1])



class readfitsobj:
	
	def __init__(self, date_dir, DAY, BAND, OBNAME):
    
		self.oname = []
		self.flist = []

		v = DAY
		pattern1 = '[' + BAND + ']*.fits'
		pattern2 = os.path.join(date_dir, v, pattern1)
		flist1 = glob.glob(pattern2)
		flist1.sort()
    
		for fl1 in flist1:
		
			hdu = fits.open(fl1)
			object_name = hdu[0].header['OBJECT']
			hdu.close()
	    
			if object_name.find(OBNAME[-4:])  != -1:
				self.oname.append(object_name)
				self.flist.append(fl1)


class quickobj:

	def __init__(self, date_dir, DAY):
    
		self.oname = []
		self.flist = []

		v = DAY
		pattern1 = '[jhk]*.fits'
		pattern2 = os.path.join(date_dir, v, 'quick', pattern1)
		flist1 = glob.glob(pattern2)
		flist1.sort()
    
		for fl1 in flist1:
        
			hdu = fits.open(fl1)
			object_name = hdu[0].header['OBJECT']
			hdu.close()
        
			self.oname.append(object_name)
			self.flist.append(fl1)
    
    

class readfitshead:
	
	def __init__(self, type1, pattern):
		
		self.flist = []
		self.oname = []
		self.band = []
		self.mjd = []
		self.ra = []
		self.dec = []
		self.offra = []
		self.offde = []
		self.comb = []
        
            
		for file in flist:
		
			hdu = fits.open(file)
			self.oname.append(hdu[0].header['OBJECT'])
			self.band.append(hdu[0].header['FILTER'])
			self.mjd.append(hdu[0].header['EXPTIME'])
			self.ra.append(hdu[0].header['RA'])
			self.dec.append(hdu[0].header['DEC'])
			self.offra.append(hdu[0].header['OFFSETRA'])
			self.offde.append(hdu[0].header['OFFSETDE'])
			self.comb.append(hdu[0].header['IMCMB001'])
			self.comb.append(hdu[0].header['IMCMB002'])
			self.comb.append(hdu[0].header['IMCMB003'])
			self.comb.append(hdu[0].header['IMCMB004'])
			self.comb.append(hdu[0].header['IMCMB005'])
			self.comb.append(hdu[0].header['IMCMB006'])
			self.comb.append(hdu[0].header['IMCMB007'])
			self.comb.append(hdu[0].header['IMCMB008'])
			self.comb.append(hdu[0].header['IMCMB009'])
			hdu.close()





def execute():
    
#copy img
#cut img
#flat
#make sky
#arith
