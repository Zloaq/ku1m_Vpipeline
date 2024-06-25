#!/Users/motomo/opt/anaconda3/bin/python3

import os
import sys
import re

from astropy.io import fits

import bottom
import dostack_auto



def comb_pset(fitslist, day, obname):
    #自分の好きなobject の名前を反映させる。2024/06/12

    header = bottom.readheader(fitslist)
    unique_band = set(header.band)
    unique_obname = set(header.object)
    uni_band = list(unique_band)
    uni_obname = list(unique_obname)


    sep_bandset = {}
    for band in uni_band:
        for obname in uni_obname:
            sep_bandset[band+obname] = [f'{band}{day}_{obname}.fits']
    for index, fitsname in enumerate(fitslist):
        band = header.band[index]
        obname = header.object[index]
        sep_bandset[band+obname].append(fitsname)
    
    for key in sep_bandset:
        bottom.combine(sep_bandset[key][1:], sep_bandset[key][0], 'average')


def comb_all(fitslist, day, obname):

    sep_band = {}
    header = bottom.readheader(fitslist)
    unique_band = set(header.band)
    uni_band = list(unique_band)

    for band in uni_band:
        sep_band[band] = [f'{band}{day}_{obname}.fits']
    for index, fitsname in enumerate(fitslist):
        band = header.band[index]
        sep_band[band].append(fitsname)

    for key in sep_band:
        bottom.combine(sep_band[key][1:], sep_band[key][0], 'average')


def match_star_coords(fitslist, refcoofile, fwhms=[3, 4, 5]):
    print('match')
    # Reference_bands.coo
    # daofind のために半値幅を選定
    # phot のための半値幅
    for fitsname in fitslist:
        data = fits.getdata(fitsname)
        nxblock = len(data[0])
        nyblock = len(data)
        starnum = {}
        for fwhm in fwhms:
            threshold = 5
            while threshold > 3:
                outputf = re.sub('.fits', f'{fwhm}.coo', fitsname)
                os.rm(outputf)
                bottom.daofind(fitsname, outputf, fwhm, threshold)
                dictlist = dostack_auto.daofind_result([outputf], [fitsname], nxblock, nyblock)
                if len(dictlist['xcen']) < 3:
                    threshold -= 0.5
                    continue
                starnum[fwhm] = len(dictlist['xcen'])
                break
            if not starnum[fwhm]:
                print(f'no star {fitsname} fwhm={fwhm}')
                break



def phot(fitslist, coordslist, fwhm):
    print('phot')

def method1():
    print('method1')

if __name__ == '__main__':
    print('mosule')