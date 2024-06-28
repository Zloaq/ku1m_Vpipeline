#!/Users/motomo/opt/anaconda3/bin/python3


import os
import sys
import re
import glob
import shutil
import subprocess
import warnings
import statistics
import numpy as np
from astropy.io import fits
from pyraf import iraf

import bottom



def flat_division(inlist):
    
    flats = {}
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for i1 in 'jhk':
        flat_path = os.path.join(script_dir, 'lib', f'{i1}*.fits')
        flats[i1] = glob.glob(flat_path)
        
    for i2 in inlist:
        print(i2, '/', flats[0],  re.sub(r'.fits', r'_fl.fits', i2))
        bottom.imarith(i2, '/', flats[i2[0]],  re.sub(r'.fits', r'_fl.fits', i2))



def method2_1(fitslist):
    # skylevel subtraction
    skylevlist = bottom.imstat(fitslist, 'midpt')

    levlist = [float(value) for value in skylevlist]

    for index, f2 in enumerate(fitslist):
        data, hdr = fits.getdata(f2, header=True)
        f3 = re.sub(r'.fits', r'_lev.fits', f2)
        data -= levlist[index]
        hdr['MEDIAN']=levlist[index]
        print(f2, '- ', levlist[index], '\t=', f3)
        fits.writeto(f3, data, hdr, overwrite=True)

def method2_2(fitslist):
    skylevlist = bottom.imstat(fitslist, 'midpt')
    levlist = [float(value) for value in skylevlist]
    levmean = statistics.mean(levlist)
    levratio = [levmean/varr for varr in levlist]

    for index, fitsname in enumerate(fitslist):
        data, hdr = fits.getdata(fitsname, header=True)
        fitsname2 = re.sub(r'.fits', r'_ylev.fits', fitsname)
        data = data * levratio[index]
        hdr['MEDIAN']=levlist[index]
        print(fitsname, '* ', levratio[index], '\t=', fitsname2)
        fits.writeto(fitsname2, data, hdr, overwrite=True)
    #2024/05/21 ああああああああああああ


def method3(flist, obnamelist):
    # make sky
    
    def self_sky(ondict):
        
        for key in ondict:
            out = '_' + key + '_skyimg.fits'
            bottom.combine(ondict[key], out, 'median')

    def off_sky(firstoff, firston, offdict):

        for band in firstoff:
            off_hdulist = [bottom.readheader(fits) for fits in firstoff[band]]
            on_hdulist  = [bottom.readheader(fits) for fits in firston[band]]
            for hdu in on_hdulist:
                base = float(hdu.mjd)
                sky_hdu = min(off_hdulist, key=lambda x: abs(float(x.mjd) - base))
                out = '_' + band + '_' + hdu.object + '_skyimg.fits'
                bottom.combine(offdict[band+'_'+hdu.object], out, 'median', 'none')
                print(out, 'was made')
        
    firstoff = {}
    firston  = {}
    offdict = {}
    ondict  = {}
    
    for index, obname in enumerate(obnamelist):
        band = flist[index][0]
        
        if 'sky' in obname:
            if band+'_'+obname not in offdict:
                offdict[band+'_'+obname] = []
                if band not in firstoff:
                    firstoff[band] = []
                    
                firstoff[band].append(flist[index])
            offdict[band+'_'+obname].append(flist[index])
            
        else:
            if band+'_'+obname not in ondict:
                ondict[band+'_'+obname] = []
                if band not in firston:
                    firston[band] = []
                    
                firston[band].append(flist[index])
            ondict[band+'_'+obname].append(flist[index])

    if len(offdict) == 0:
        self_sky(ondict)
    else:
        off_sky(firstoff, firston, offdict)
        
    return [item for sublist in ondict.values() for item in sublist]



def method4(flist, obnamelist):
    # sky subtraction

    for i1 in range(len(obnamelist)):
        band = flist[i1][0]
        sky = '_' + band + '_' + obnamelist[i1] + '_skyimg.fits'
        out2 = re.sub(r'.fits', r'_sky.fits', flist[i1])
        bottom.imarith(flist[i1], '-', sky, out2)
        print(flist[i1],'-', sky, '=' ,out2)



def execute(day, band, obname):

    if param.quicklook == 1:
        
        print('start quicklook')
        param.cut = 1
        param.mksky = 1
        param.skysub = 1
        print("quicklook is not available.")
        sys.exit()
        
    elif param.quicklook == 0:
        
        print('serch fits of',argvs[2])
        method1()
        os.makedirs(pattern2, exist_ok=True)
    
    
    if param.cut == 1:
        # cut black frame

        print('cut fits')
        bottom.cut_copy(fitslist, '[9:328,9:264]')
        
    iraf.chdir(pattern2)


    if param.flatdiv == 1:
        flat_division()
    
    
    if param.skylevsub == 1:
        # skylevel

        print('subtraction of sky level')
        fitslist1 = method0()
        method2_1(fitslist1)
        
    
    if param.mksky == 1:
        # make sky

        print('make sky')
        fitslist1 = method0()
        method3(fitslist1)
        
    if param.skysub == 1:
        fitslist1 = method0()
        method4(fitslist1)
    
    
    if param.rm_cut == 1:
        subprocess.run('rm *_cut.fits', shell=True)
    if param.rm_st == 1:
        subprocess.run('rm *_st.fits', shell=True)
    if param.rm_sky	== 1:
        subprocess.run('rm *_sky.fits', shell=True)
    
    write_log()


    if param.dostack_a == 1:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_dir, 'dostack_auto.py')
        subprocess.run([script_path, argvs[1], argvs[2], argvs[3]])
        
        
    if param.phot == 1:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_dir, 'aperture_phot.py')
        subprocess.run([script_path, argvs[1], argvs[2], argvs[3]])



if __name__ == "__main__":

    argvs = sys.argv
    argc = len(argvs)
    obnamelist = []
    fitslist = []
    fitsid = []
    fitspro = []

    param = readparam()
    pattern1 = os.path.join(param.date_dir, argvs[1])
    pattern2 = os.path.join(param.date_dir, argvs[1], argvs[2])
    
    iraf.chdir(pattern1)
    
    list1 = glob.glob('*-????.fits')
    list1.sort()
    ist1 = readheader(list1)
	
    if argc ==  4:
        execute(argvs[1], argvs[3], argvs[2])
        
    elif argc == 2:#yetyetyet
        execute(argvs[1])
        
    else:
        print('iraf_script.py [DAY] [obname] [band]')
        print('usage: iraf_script.py 230103 ZTF23aaaatwl jhk')
        sys.exit()

