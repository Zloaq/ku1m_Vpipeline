#!/Users/motomo/opt/anaconda3/bin/python3

import os
import re
import sys
import math
import scipy as sp
import numpy as np
from astropy.io import fits
from scipy.stats import mode
from pyraf import iraf



class readheader():
    
    def __init__(self, inlist):
        
        bandset = {'g':'g', 'i':'i', 'J':'j', 'H':'h', 'Ks':'k'}

        self.filename = []
        self.object = []
        self.exptime = []
        self.band = []
        self.filter = []
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
        
        for f1 in inlist:
            
            hdu = fits.open(f1)
            self.filename.append(f1)
            self.object.append(hdu[0].header['OBJECT'])
            self.exptime.append(hdu[0].header['EXPTIME'])
            self.band.append(bandset[hdu[0].header['FILTER']])
            self.filter.append(hdu[0].header['FILTER'])
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



def setparam(band='j'):
	
    scale = {
         'g': 0.692, 'i': 0.692, 'j': 0.692, 'h': 0.692, 'k': 0.692
    }
    readnoise = {
        'g': 11.97, 'i': 11.97, 'j': 11.97, 'h': 11.97, 'k': 11.97
    }
    epadu = {
        'g': 11.97, 'i': 11.97, 'j': 10, 'h': 11.7, 'k': 9.5
    }

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
    iraf.apphot.datapars.scale = scale[band]
    iraf.apphot.datapars.fwhmpsf = 3
    iraf.apphot.datapars.datamax = 20000
    iraf.apphot.datapars.datamin = 0
    iraf.apphot.datapars.readnoise = readnoise[band]
    iraf.apphot.datapars.epadu = epadu[band]
    iraf.apphot.centerpars.calgorithm = 'centroid'
    iraf.apphot.fitskypars.salgorithm = 'mode'
    iraf.apphot.fitskypars.dannulus = 10
    iraf.apphot.findpars.threshold = 5
    iraf.apphot.findpars.sharphi = 0.8



def cut_copy(inlist, gicrange, jhkcrange):
    
    with open('_inlist', 'w') as fin, open('_outlist', 'w') as fout:
        cutrange = {
            'g':gicrange, 'i':gicrange, 'j':jhkcrange, 'h':jhkcrange, 'k':jhkcrange
                    }
        for f1 in inlist:
            crange = cutrange[inlist[0]]
            f2 = re.sub(r'.fits', '_cut.fits', f1)
            f3 = f1+crange
            fin.write(f3 + '\n')
            fout.write(f2 + '\n')
            print('imcopy', f1, crange, f2)
    
    iraf.imcopy(input='@_inlist', output='@_outlist', Stdout=1)



def imstat(inlist, method):
    
    with open('_inlist', 'w')as f1:
        for i1 in inlist:
            f1.write(i1+'\n')
    
    iraf.imstatistics.unlearn()
    iraf.imstatistics.fields = method
    iraf.imstatistics.nclip = 5
    iraf.imstatistics.format = 'No'
    return iraf.imstatistics('@_inlist', Stdout=1)


def skystat(fitsname, method, nclip=3, clip_threshold=3):
    
    data = fits.getdata(fitsname)
    flat_data = data.flatten()
    threshold = np.percentile(flat_data, 100)
    backgroud_data = flat_data[flat_data <= threshold]
    mode_value, count = mode(backgroud_data)
    median_value = np.median(backgroud_data)

    if method == 'mode':
        return mode_value
    
    if method == 'mode,stddev':
        basis_value = mode_value
    elif method == 'median,stddev':
        basis_value = median_value
    else:
        basis_value = median_value

    for num in range(nclip):
        squared_diffs = (backgroud_data - basis_value) ** 2
        mean_squared_diff = np.mean(squared_diffs)
        rms_from_mode = np.sqrt(mean_squared_diff)
        lower_bound = basis_value - clip_threshold * rms_from_mode
        upper_bound = basis_value + clip_threshold * rms_from_mode
        backgroud_data = backgroud_data[(backgroud_data >= lower_bound) & (backgroud_data <= upper_bound)]
    skyrms = rms_from_mode
    print('うおおおおおおおおおおおお', basis_value, skyrms)

    if method == 'stddev':
        return skyrms
    elif method == 'median':
        return median_value
    elif method == 'mode,stddev':
        return basis_value, skyrms
    elif method == 'median,stddev':
        return basis_value, skyrms


def combine(inlist, outname, method, reject='sigclip'):
	
    with open('_inlist', 'w')as f1:
        for i1 in inlist:
            f1.write(i1+'\n')
        
    iraf.imcombine.unlearn()
    iraf.imcombine.combine = method
    iraf.imcombine.reject = reject
    iraf.imcombine.zero = 'none'
    iraf.imcombine.mclip = 'yes'
    iraf.imcombine.lsigma = 3.0
    iraf.imcombine.hsigma = 3.0
    
    iraf.imcombine('@_inlist', outname, Stdout=1)



def imarith(infile1, op, infile2,  outfile):

    if op != '+' and op != '-' and op != '/':
            print("bad operator, op")
            return 1

    data1, hdr = fits.getdata(infile1, header=True)
    data2 = fits.getdata(infile2)

    if op == '+':
            data3 = data1 + data2
    elif op == '-':
            data3 = data1 - data2
    elif op == '/':
            data3 = data1 / data2

    fits.writeto(outfile, data3, hdr, overwrite=True)



def daofind(inputf, outputf, fwhm, threshold=5 ):
    
    iraf.noao()
    iraf.digiphot(Stdout=1)
    iraf.apphot(Stdout=1)

    mode, stddev = skystat(inputf, 'mode,stddev')

    iraf.apphot.datapars.scale = 0.692
    iraf.apphot.datapars.fwhmpsf = fwhm
    iraf.apphot.datapars.sigma = stddev
    iraf.apphot.datapars.datamin = float(mode)-float(stddev)
    iraf.apphot.datapars.datamax = 15000
    iraf.apphot.findpars.threshold = threshold
    iraf.daofind.unlearn()
    iraf.daofind.output = 'temp1.coo'
    iraf.daofind.interac = 'no'
    iraf.daofind.verify = 'no'
    iraf.daofind.mode = 'h'
    iraf.daofind(inputf)

    iraf.apphot.center.unlearn()
    iraf.apphot.centerpars.calgorithm = 'centroid'
    iraf.apphot.center.coords = 'temp1.coo'
    iraf.apphot.center.output = 'temp2.coo'
    iraf.apphot.center.interactive = 'no'
    iraf.apphot.center.verify = 'no'
    iraf.apphot.center.update = 'yes'
    iraf.apphot.center.verbose = 'no'
    iraf.apphot.center.mode = 'h'
    iraf.apphot.center(inputf)

    if not os.path.exists('temp2.coo'):
        print(f'temp2.coo is not exists')
        file2 = open(outputf, 'w')
        file2.close()
        os.remove('temp1.coo')
        return

    with open('temp2.coo', 'r') as file1:
        lines = file1.readlines()
    filtered_lines = [line for line in lines if not "temp1.coo" in line]
    with open(outputf, 'w') as file2:
        file2.writelines(filtered_lines)

    os.remove('temp1.coo')
    os.remove('temp2.coo')



def psfmeasure(inputf):
	
    iraf.noao()
    iraf.obsutil(Stdout=1)
    iraf.psfmeasure.unlearn()
    iraf.psfmeasure.coords = 'center'
    iraf.psfmeasure.scale = 0.69
    iraf.psfmeasure.display = 'No'
    iraf.psfmeasure.logfile = ''
    iraf.psfmeasure.mode = 'h'
    
    return iraf.psfmeasure(inputf, Stdout=1)



def imshift(inputf, outputf, xshift, yshift):
	
	iraf.imshift.mode = 'h'
	iraf.imshift(inputf, outputf, xshift, yshift)


def xyxymatch(inputf, referencef, outputf):
    
	iraf.xyxymatch.unlearn()
	iraf.xyxymatch.input = inputf
	iraf.xyxymatch.reference = referencef
	iraf.xyxymatch.output = outputf
	iraf.xyxymatch.tolerance = 3.0
	iraf.xyxymatch.matching = 'triangles'
	iraf.xyxymatch.interactive = 'No'
	iraf.xyxymatch.verbose = 'Yes'
	iraf.xyxymatch.mode = 'h'
	iraf.xyxymatch(inputf)


def lnum_match(inputf, referancef, outputf): # 変な画像があることは考慮してない。

    def read_coofile(infile):
        with open(infile, 'r') as file:
            lines = [line.strip().split() for line in file.readlines() if not line.startswith('#')]
        coords = [[float(line[0]), float(line[1]), float(line[2])] for line in lines] # coox, cooy, mag
        return coords

    def calculate_segments(coords):
        segments = []
        points = len(coords)
        if points <= 1:
            return segments
        for p1 in range(points-1):
            x1 = coords[p1][0]
            y1 = coords[p1][1]
            for p2 in range(p1+1, points):
                x2 = coords[p2][0]
                y2 = coords[p2][1]

                length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
                angle = math.atan2(y2-y1, x2-x1)

                if y1 < y2:
                    segments.append([(x1, y1), (x2, y2), length, angle])
                else:
                    segments.append([(x2, y2), (x1, y1), length, angle])

        return segments


    def match_segments(inlist, reflist, length_tolerance=1e-2, angle_tolerance=1e-2):
        print('match')
        match_coo = []
        for seg1 in reflist:
            rcoo1 = seg1[0]
            rcoo2 = seg1[1]
            rlength = seg1[2]
            rangle = seg1[3]
            for seg2 in inlist:
                icoo1 = seg2[0]
                icoo2 = seg2[1]
                ilength = seg2[2]
                iangle = seg2[3]

                length_diff = abs(rlength - ilength)
                angle_diff = abs(rangle - iangle)
                if length_diff <= length_tolerance and angle_diff <= angle_tolerance:
                    match_coo.append((rcoo1, icoo1))
                    match_coo.append((rcoo2, icoo2))
                    break
        return match_coo
    

    def one_point_match(incoords, refcoords):#一番明るさの近いやつでいいや
        print('one point match')
        difflist = []
        for index1 in range(len(refcoords)):
            x1 = refcoords[index1][0]
            y1 = refcoords[index1][1]
            m1 = refcoords[index1][2]
            for index2 in range(len(incoords)):
                x2 = incoords[index2][0]
                y2 = incoords[index2][1]
                m2 = incoords[index2][2]

                magdiff = abs(m1 - m2)
                difflist.append([(x1, y1), (x2, y2), magdiff])

        match_coo = min(difflist, key=lambda x:x[2])

        return [match_coo]

    #途中ううううううう 2024/05/16
    refcoords = read_coofile(referancef)
    inpcoords = read_coofile(inputf)

    if len(refcoords) == 0 or len(inpcoords) == 0:
        raise Exception('座標がありません')

    refseg = calculate_segments(refcoords)
    inpseg = calculate_segments(inpcoords)

    if len(refseg) == 0 or len(inpseg) == 0:
        match_coo = one_point_match(inpcoords, refcoords)
    else:
        match_coo = match_segments(inpseg, refseg)

    # write .match
    with open(outputf, 'w') as f1:
        f1.write(
            f'# Input: {inputf}\n'
            f'# Reference: {referancef}\n'
            f'# Column definitions\n'
            f'#    Column 1: X reference coordinate\n'
            f'#    Column 2: Y reference coordinate\n'
            f'#    Column 3: X input coordinate\n'
            f'#    Column 4: Y input coordinate\n\n'
        )

        print('testtttttttttttttttt',match_coo)
        for coo in match_coo:
            print('testttttttttt2', coo)
            f1.write(
                f'   '
                '{:<7}'.format(coo[0][0])+'   '
                '{:<7}'.format(coo[0][1])+'   '
                '{:<7}'.format(coo[1][0])+'   '
                '{:<7}'.format(coo[1][1])+'\n'
                )

def geomap(inputf, xmaxv, ymaxv, databasef, fitmetry):
    
    iraf.geomap.unlearn()
    iraf.geomap.input = inputf
    iraf.geomap.database = databasef
    iraf.geomap.xmin = 1
    iraf.geomap.xmax = xmaxv
    iraf.geomap.ymin = 1
    iraf.geomap.ymax = ymaxv
    iraf.geomap.fitgeometry = fitmetry
    iraf.geomap.interactive = 'No'
    iraf.geomap.mode = 'h'
    iraf.geomap(inputf, databasef, Stdout=1)



def geotran_param(xmaxv, ymaxv):
    
    iraf.geotran.unlearn()
    
    iraf.geotran.database = ''
    iraf.geotran.transforms = ''
    iraf.geotran.geometry = 'linear'
    iraf.geotran.xmin = 1
    iraf.geotran.xmax = xmaxv
    iraf.geotran.ymin = 1
    iraf.geotran.ymax = ymaxv
    iraf.geotran.xscale = 1.
    iraf.geotran.yscale = 1.
    iraf.geotran.interpo = 'drizzle'
    iraf.geotran.boundary = 'nearest'
    iraf.geotran.constant = -100000
    iraf.geotran.fluxconserve = 'yes'
    iraf.geotran.nxblock = xmaxv
    iraf.geotran.nyblock = ymaxv
    iraf.geotran.verbose = 'no'
    iraf.geotran.mode = 'h'


def phot(fitsname, fwhm, sigma):

    iraf.imstatistics.unlearn()
    iraf.imstatistics.fields = 'mode,stddev'
    iraf.imstatistics.nclip = 5
    iraf.imstatistics.format = 'No'
    mode, stddev = iraf.imstatistics(fitsname, Stdout=1)
    
    iraf.apphot.datapars.fwhmpsf = fwhm
    iraf.apphot.datapars.sigma = sigma
    iraf.apphot.datapars.datamin = float(mode)-2*float(stddev)
    iraf.apphot.datapars.datamax = 15000
    iraf.apphot.datapars.readnoise = 1
    iraf.apphot.datapars.epadu = 1
    #まだだだだだだ

    iraf.apphot.fitskypars.salgorithm = 'mode'
    iraf.apphot.fitskypars.annulus = 2.5*float(fwhm)
    iraf.apphot.fitskypars.dannulus = 10
    iraf.apphot.photpars.aperttures = fwhm
    iraf.apphot.photpars.zmag = 0

    iraf.apphot.phot.interactive = 'no'
    iraf.apphot.phot.cache = 'no'
    iraf.apphot.phot.verify = 'no'
    iraf.apphot.phot.update = 'yes'
    iraf.apphot.phot.verbose = 'no'
    iraf.apphot.phot.mode = 'h'




if __name__ == '__main__':
	
	print('bottom.py is a module.')
