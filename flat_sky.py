#!/opt/anaconda3/envs/p11/bin/python3


import os
import sys
import re
import glob
import subprocess
import statistics
import numpy as np
from tqdm import tqdm
from astropy.io import fits
import matplotlib.pyplot as plt
from pyraf import iraf

import bottom



def flat_division(bands, inlist):
    
    flats = {}
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for i1 in bands:
        flat_path = os.path.join(script_dir, 'lib', f'{i1}*.fits')
        flats[i1] = glob.glob(flat_path)
        if not flats[i1]:
            return 1
        
    for i2 in tqdm(inlist, desc='{:<13}'.format('flat div')):
        bottom.imarith(i2, '/', flats[i2[0]][0],  re.sub(r'.fits', r'_fl.fits', i2))

    return 0

def method2_0(fitslist0):
    # skylevel subtraction
    for band, fitslist in fitslist0.items():
        outliers_index = []
        levlist = [bottom.skystat(fitsname, 'median') for fitsname in tqdm(fitslist, desc='{:<13}'.format(f'calc {band} skylev'))]
        if band in ['g', 'i']:
            # 標準偏差を計算
            levlist_np = np.array(levlist)
            median = np.median(levlist_np)
            std_dev = np.sqrt(np.mean((levlist_np - median) ** 2))
            lower_bound = median - 3 * std_dev
            upper_bound = median + 3 * std_dev
            # 外れ値のインデックスを取得
            outliers_index = [index for index, value in enumerate(levlist) if value < lower_bound or value > upper_bound]
        for index in outliers_index:
            print(f'{fitslist[index]} was rejected as outside 3 sigma skylevel')
        for index, f2 in enumerate(tqdm(fitslist, desc='{:<13}'.format(f'sub  {band} skylev'))):
            if index in outliers_index:
                continue  # 外れ値の場合はスキップ
            data, hdr = fits.getdata(f2, header=True)
            f3 = re.sub(r'.fits', r'_lev.fits', f2)
            data = data.astype(float)
            data -= levlist[index]
            hdr['SKYCOUNT'] = levlist[index]
            fits.writeto(f3, data, hdr, overwrite=True)



def method2_1(fitslist):
    # skylevel subtraction
    levlist = [bottom.skystat(fitsname, 'median') for fitsname in tqdm(fitslist, desc='{:<13}'.format('calc skylevel'))]
    for index, f2 in enumerate(tqdm(fitslist, desc='{:<13}'.format('sub skylevel'))):
        data, hdr = fits.getdata(f2, header=True)
        f3 = re.sub(r'.fits', r'_lev.fits', f2)
        #print(f'{f2} {levlist}')
        data = data.astype(float)
        data -= levlist[index]
        hdr['SKYCOUNT']=levlist[index]
        #print(f2, '- ', levlist[index], '\t=', f3)
        fits.writeto(f3, data, hdr, overwrite=True)

def method2_2(fitslist):
    levlist = [bottom.skystat(fitsname, 'median') for fitsname in fitslist]
    levmean = statistics.mean(levlist)
    levratio = [levmean/varr for varr in levlist]

    for index, fitsname in enumerate(tqdm(fitslist, desc='{:<13}'.format('div skylevel'))):
        data, hdr = fits.getdata(fitsname, header=True)
        fitsname2 = re.sub(r'.fits', r'_ylev.fits', fitsname)
        data = data * levratio[index]
        hdr['SKYCOUNT']=levlist[index]
        #print(fitsname, '* ', levratio[index], '\t=', fitsname2)
        fits.writeto(fitsname2, data, hdr, overwrite=True)
    #2024/05/21 ああああああああああああ


def method3(flist, obnamelist):
    # make sky
    
    def self_sky(ondict):
        
        for key in tqdm(ondict, desc='{:<13}'.format('make selfsky')):
            out = '_' + key + '_skyimg.fits'
            bottom.combine(ondict[key], out, 'median')

    def off_sky(firstoff, firston, offdict):

        for band in tqdm(firstoff, desc='{:<13}'.format('make offsky')):
            #print(f'\nband {band}')
            #print(f'\nfirstoff{firstoff}')
            #print(f'\nfirston{firston}')
            off_hduobj = bottom.readheader(firstoff[band])
            on_hduobj  = bottom.readheader(firston[band])
            for idxon, mjd in enumerate(on_hduobj.mjd):
                idxoff = min(range(len(off_hduobj.mjd)), key=lambda i: abs(float(off_hduobj.mjd[i]) - mjd))
                out = '_' + band + '_' + on_hduobj.object[idxon] + '_skyimg.fits'
                #print(f'check {band+"_"+off_hduobj.object[idxoff]}')
                #print(offdict)
                bottom.combine(offdict[band+'_'+off_hduobj.object[idxoff]], out, 'median', 'none')
                #print(out, 'was made')
        
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
        skytype = 'self'
    else:
        off_sky(firstoff, firston, offdict)
        skytype = 'off'
        
    return [item for sublist in ondict.values() for item in sublist], skytype



def method4(flist, obnamelist, skytype):
    # sky subtraction

    for i1 in tqdm(range(len(obnamelist)), desc='{:<13}'.format('sky_subtract')):
        band = flist[i1][0]
        sky = '_' + band + '_' + obnamelist[i1] + '_skyimg.fits'
        out2 = re.sub(r'.fits', f'_{skytype}sky.fits', flist[i1])
        bottom.imarith(flist[i1], '-', sky, out2)
        #print(flist[i1],'-', sky, '=' ,out2)



if __name__ == "__main__":

	print('flat_sky.py is a module.')