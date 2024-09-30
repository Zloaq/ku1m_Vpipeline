#!/opt/anaconda3/envs/p11/bin/python3

import os
import re
import sys
import math
import subprocess
import statistics
import numpy as np
import pandas as pd
from scipy import ndimage
from tqdm import tqdm
from astropy.io import fits
from pyraf import iraf

from itertools import combinations
from scipy.spatial import KDTree
from scipy.spatial.transform import Rotation as R
from sklearn.neighbors import NearestNeighbors

import bottom

import matplotlib.pyplot as plt

def starfind_center3(fitslist, pixscale, satcount, searchrange=[3.0, 5.0, 0.2], minstarnum=0, maxstarnum=100):
    
    def squareness(region_slice):
        width = region_slice[1].stop - region_slice[1].start
        height = region_slice[0].stop - region_slice[0].start
        squareness = abs(width - height)
        if width < 3 or height < 3:
            squareness = 20
        return squareness

    def filling_rate(label, region_slice, labeled_image):
        width = region_slice[1].stop - region_slice[1].start
        height = region_slice[0].stop - region_slice[0].start
        #print(f'region_slice\n{region_slice}')
        #sys.exit()
        region = labeled_image[region_slice]
        area = np.sum(region == label)
        expected_area = width * height
        filling_ratio = area / expected_area
        return filling_ratio

    def filter_data(data, threshold, rms):
        return np.array([[d if d > threshold * rms else 0 for d in r] for r in data])
    
    def binarize_data(data, threshold, rms):
        return np.array([[1 if d > threshold * rms else 0 for d in r] for r in data])
    
    def filter_saturate(labeled_image, filtered_data, objects, header):
        #div には対応してない
        skycount = float(header['SKYCOUNT'] or 0)
        saturation_value = float(satcount) - skycount
        saturated_mask = filtered_data >= saturation_value
        saturated_labels = np.unique(labeled_image[saturated_mask])
        
        templabels = np.arange(1, len(objects) + 1)
        mask = ~np.isin(templabels, saturated_labels)
        filtered_labels = templabels[mask]
        filtered_objects = np.array(objects)[mask]
        filtered_objects = [tuple(inner_list) for inner_list in filtered_objects.tolist()]
        #print(f'filtered_labels\n{filtered_labels.tolist()}')
        #print(f'filtered_objects\n{tuple(filtered_objects.tolist())}')
        #sys.exit()
        
        return filtered_labels.tolist(), filtered_objects

    def detect_round_clusters(filtered_labels, filtered_objects, labeled_image, square=3, fillrate=0.5):
        squareness_values = np.array([squareness(region_slice) for region_slice in filtered_objects])
        filling_rates = np.array([filling_rate(filtered_labels[i], region_slice, labeled_image) 
                                for i, region_slice in enumerate(filtered_objects)])

        mask = (squareness_values < square) & (filling_rates > fillrate)
        
        round_clusters = np.array(filtered_labels)[mask]
        slice_list = np.array([[region_slice[0], region_slice[1]] for region_slice in filtered_objects])[mask]

        return round_clusters, slice_list
    
    def clustar_center(slices):
        centers = []
        for sl in slices:
            y_slice, x_slice = sl
            if x_slice.stop - x_slice.start < 3 or y_slice.stop - y_slice.start < 3:
                continue
            center_row = (y_slice.start + y_slice.stop - 1) / 2
            center_col = (x_slice.start + x_slice.stop - 1) / 2
            centers.append([center_row, center_col])
        return centers

    def write_to_txt(centers, filename):
        with open(filename, 'w') as f1:
            for center in centers:
                f1.write(f'{center[1]}  {center[0]}\n')

    def chose_unique_coords(coordsfile):
        input_file_path = coordsfile
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        header_lines = []
        data_lines = []
        for line in lines:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                data_lines.append(line)
        """
        data = []
        for line in data_lines:
            parts = line.strip().split()
            if len(parts) == 8:
                data.append(parts)
        """
        data = [line.strip().split() for line in data_lines if len(line.strip().split()) == 8]


        df = pd.DataFrame(data, columns=['XCENTER', 'YCENTER', 'XSHIFT', 'YSHIFT', 'XERR', 'YERR', 'CIER', 'CERROR'])
        df = df.astype({'XCENTER': 'float', 'YCENTER': 'float', 'XSHIFT': 'float', 'YSHIFT': 'float', 'XERR': 'float', 'YERR': 'float', 'CIER': 'int', 'CERROR': 'str'})

        df['X_INT'] = df['XCENTER'].astype(int)
        df['Y_INT'] = df['YCENTER'].astype(int)
        df['COUNT'] = df.groupby(['XCENTER', 'YCENTER'])['XCENTER'].transform('count')
        idx = df.groupby(['X_INT', 'Y_INT'])['COUNT'].idxmax()
        result = df.loc[idx, ['XCENTER', 'YCENTER', 'XSHIFT', 'YSHIFT', 'XERR', 'YERR', 'CIER', 'CERROR']].drop_duplicates(subset=['XCENTER', 'YCENTER'])

        starnum = 0
        with open(coordsfile, 'w') as file:
            for line in header_lines:
                file.write(line)
            result_reset = result.reset_index(drop=True)
            for index, row in result_reset.iterrows():
                line = f"{row['XCENTER']:<14.3f} {row['YCENTER']:<11.3f} {row['XSHIFT']:<8.3f} {row['YSHIFT']:<8.3f} {row['XERR']:<8.3f} {row['YERR']:<15.3f} {row['CIER']:<5d} {row['CERROR']:<9s}\n"
                file.write(line)
                starnum = index + 1

        return starnum, coordsfile

    
    coordsfilelist = []
    starnumlist = []
    threshold_lside = []
    
    #print(f'{band} band threshold range = {searchrange}')

    #for index, filename in enumerate(tqdm(fitslist, desc=f'{fitslist[0][0]} band starfind')):
    iterate = 0
    searchrange0 = searchrange
    with tqdm(total=len(fitslist), desc=f'{fitslist[0][0]} band starfind') as pbar:
        while iterate <= len(fitslist) - 1:
            filename = fitslist[iterate]
            data = fits.getdata(filename)
            header = fits.getheader(filename)
            offra_pix = int(float(header['OFFSETRA'])/pixscale)
            offde_pix = int(float(header['OFFSETDE'])/pixscale)
            if offra_pix > 0:
                data[:, :offra_pix] = 0
            elif offra_pix < 0:
                data[:, offra_pix:] = 0
            if offde_pix > 0:
                data[-offde_pix:, :] = 0
            elif offde_pix < 0:
                data[:-offde_pix, :] = 0
            
            rms = bottom.skystat(filename, 'stddev')

            roopnum = [0, 0]
            while searchrange0[0] >= 2.0:
                center_list = []
                for threshold in np.arange(searchrange0[0], searchrange0[1], searchrange0[2]):
                    binarized_data = binarize_data(data, threshold, rms)
                    filtered_data = filter_data(data, threshold, rms)
                    labeled_image, _ = ndimage.label(binarized_data)
                    objects = ndimage.find_objects(labeled_image)
                    filtered_labels, filtered_objects = filter_saturate(labeled_image, filtered_data, objects, header)

                    _, slice_list = detect_round_clusters(filtered_labels, filtered_objects, labeled_image)

                    centers = clustar_center(slice_list)
                    center_list.extend(centers)
                write_to_txt(center_list, 'temp.coo')
                coordsfile = re.sub('.fits', '.coo', filename)
                #これで時間食ってる
                bottom.center(filename, coordsfile, 'temp.coo')
                starnum, file = chose_unique_coords(coordsfile)
                #print(f'{filename}, {searchrange0}')

                if roopnum[0] > 0 and roopnum[1] > 0:
                    if searchrange0[0] > 2.0:
                        searchrange0[0] -= 1
                        searchrange0[1] -= 1
                    break

                elif (starnum < minstarnum) & (searchrange0[0] > 2.0):
                    #print(f'retry starfind')
                    searchrange0[0] -= 1
                    searchrange0[1] -= 1
                    roopnum[0] += 1
                    continue
                elif starnum > maxstarnum:
                    searchrange0[0] += 1
                    searchrange0[1] += 1
                    roopnum[1] += 1
                    continue
                else:
                    break

            coordsfilelist.append(file)
            starnumlist.append(starnum)
            iterate += 1
            threshold_lside.append(searchrange0[0])
            pbar.update(1)

    return starnumlist, coordsfilelist, threshold_lside


def triangle_match(inputf, referencef, outputf, match_threshold=0.05):
    def read_coofile(infile):
        with open(infile, 'r') as file:
            flines = file.readlines()
        lines = [line.strip().split() for line in flines if not line.startswith('#')]
        coords = np.array([[float(line[0]), float(line[1])] for line in lines])  # coox, cooy,
        return coords

    def compute_triangle_descriptors(coords):
        triangles = list(combinations(coords, 3))
        descriptors = []
        for tri in triangles:
            tri = np.array(tri)
            # 辺の長さを計算
            sides = [np.linalg.norm(tri[i] - tri[j]) for i, j in combinations(range(3), 2)]
            sides = np.sort(sides)
            
            #sides = sides / sides[-1]  # 最大の辺の長さで割る

            a, b, c = sides
            angles = [
                np.arccos((b**2 + c**2 - a**2) / (2 * b * c)),
                np.arccos((a**2 + c**2 - b**2) / (2 * a * c)),
                np.arccos((a**2 + b**2 - c**2) / (2 * a * b))
            ]
            angles = np.sort(angles)
            descriptor = np.concatenate([sides, angles])
            descriptors.append((descriptor, tri))
        return descriptors

    def match_triangle():
        
        coords_input = read_coofile(inputf)
        coords_ref = read_coofile(referencef)

        descriptors_input = compute_triangle_descriptors(coords_input)
        descriptors_ref = compute_triangle_descriptors(coords_ref)
        #print(f'入力画像の三角形の数: {len(descriptors_input)}')
        #print(f'参照画像の三角形の数: {len(descriptors_ref)}')

        descs_input = np.array([desc[0] for desc in descriptors_input])
        descs_ref = np.array([desc[0] for desc in descriptors_ref])

        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(descs_ref)
        distances, indices = nbrs.kneighbors(descs_input)

        threshold = match_threshold
        good_matches = distances[:, 0] < threshold

        matched_triangles_input = [descriptors_input[i][1] for i in range(len(descriptors_input)) if good_matches[i]]
        matched_triangles_ref = [descriptors_ref[indices[i][0]][1] for i in range(len(descriptors_input)) if good_matches[i]]

        #print(f'入力画像の三角形のマッチした数: {len(matched_triangles_input)}')
        #print(f'参照画像の三角形のマッチした数: {len(matched_triangles_ref)}')

        src_points = []
        dst_points = []
        for tri_input, tri_ref in zip(matched_triangles_input, matched_triangles_ref):
            for pt_input, pt_ref in zip(tri_input, tri_ref):
                src_points.append(pt_input)
                dst_points.append(pt_ref)
        src_points = np.array(src_points)
        dst_points = np.array(dst_points)

        src_points_unique, indices = np.unique(src_points, axis=0, return_index=True)
        dst_points_unique = dst_points[indices]

        matched_pair = [(dst_points_unique[i], src_points_unique[i]) for i in range(len(dst_points_unique))]
        matched_pair = np.array(matched_pair)
        
        #print(f'src\n{src_points_unique}')
        #print(f'dst\n{dst_points_unique}')
        #print(f'matched\n{matched_pair}')

        return matched_pair
    
    matched_coo = match_triangle()
    # write .match

    if matched_coo.size == 0:
        return None

    with open(outputf, 'w') as f1:
        f1.write(
            f'# Input: {inputf}\n'
            f'# Reference: {referencef}\n'
            f'# Column definitions\n'
            f'#    Column 1: X reference coordinate\n'
            f'#    Column 2: Y reference coordinate\n'
            f'#    Column 3: X input coordinate\n'
            f'#    Column 4: Y input coordinate\n\n'
        )
        
        #print(f'testtttttttttttttttt, {matched_coo}, {matched_coo.shape}')
        for coo_varr in matched_coo:
            #print('testttttttttt2', coo)
            #print(coo_varr[0][0])
            #print(coo_varr[0][0].shape)
            f1.write(
                f'   '
                '{:<7}'.format(coo_varr[0][0].item())+'   '
                '{:<7}'.format(coo_varr[0][1].item())+'   '
                '{:<7}'.format(coo_varr[1][0].item())+'   '
                '{:<7}'.format(coo_varr[1][1].item())+'\n'
                )
    return outputf



def match_checker(checklist):
    #print('check すっぞ')
    linenumlist = []
    for file in checklist:
        with open(file, 'r') as f1:
            lines = f1.readlines()
        linenum = 0
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            linenum += 1
        linenumlist.append(linenum)
    matched_star = min(linenumlist) if linenumlist else 0

    return matched_star


def do_trimatch(optcoolist, infcoolist):
    for varr in optcoolist:
        for filename in optcoolist[varr][1:]:
            outf = re.sub(r'.coo', r'.match', filename)
            triangle_match(filename, optcoolist[varr][0], outf)


def geotparam(param, file_list, base_rotate):

    param_list = []

    for filename in file_list:
        geotp = {}
        geotp['fitsid']=filename[1:-4]
        tempfits = f'{filename[0:-4]}.fits'
        hdu = fits.open(tempfits)
        move_rotate = float(hdu[0].header['OFFSETRO']) or 0
        rotate_diff = abs(base_rotate - move_rotate)
        rotate1 = 360 - rotate_diff
        rotate2 = rotate_diff


        with open(filename, 'r') as f1:
            lines = f1.readlines()
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            line_list = line.split()
            if len(line_list) == 1:
                continue
            if 'xrefmean' in line:
                geotp['xrefmean']=float(line_list[1])
            elif 'yrefmean' in line:
                geotp['yrefmean']=float(line_list[1])
            elif 'xmean' in line:
                geotp['xmean']=float(line_list[1])
            elif 'ymean' in line:
                geotp['ymean']=float(line_list[1])
            elif 'geometry' in line:
                geotp['geometry']=line_list[1]
            elif 'xshift' in line:
                geotp['xshift']=float(line_list[1])
            elif 'yshift' in line:
                geotp['yshift']=float(line_list[1])
        
            elif 'xrotation' in line:
                if (abs(float(line_list[1]) - rotate1) < 5) or (abs(float(line_list[1]) - rotate2) < 5):
                    geotp['xrotation']=float(line_list[1])
            elif 'yrotation' in line:
                if (abs(float(line_list[1]) - rotate1) < 5) or (abs(float(line_list[1]) - rotate2) < 5):
                    geotp['yrotation']=float(line_list[1])
        
        if len(geotp) != 10:
            continue

        if geotp['xrotation'] != geotp['yrotation']:
            continue
        
        if len(geotp) == 10:
            param_list.append(geotp)
                    
    return param_list


def do_starfind(fitslist, param, optkey, infrakey):

    optstarlist = {}
    optcoolist = {}
    infstarlist = {}
    infcoolist = {}
    opt_l_threshold = {}
    inf_l_threshold = {}

    def iterate_part(fitslist0, param, h_threshold=10, l_threshold=9, interval=0.2):
        band = fitslist0[0][0]
        pixscale = {
        'g':param.pixscale_g, 'i':param.pixscale_i,
        'j':param.pixscale_j, 'h':param.pixscale_h, 'k':param.pixscale_k
        }
        satcount = {
        'g':param.g_satcount, 'i':param.i_satcount,
        'j':param.j_satcount, 'h':param.h_satcount, 'k':param.k_satcount
        }
        maxstar = {
        'g':param.g_maxstarnum, 'i':param.i_maxstarnum,
        'j':param.j_maxstarnum, 'h':param.h_maxstarnum, 'k':param.k_maxstarnum
        }
        minstar = {
        'g':param.g_minstarnum, 'i':param.i_minstarnum,
        'j':param.j_minstarnum, 'h':param.h_minstarnum, 'k':param.k_minstarnum
        }
        minstarnum = 20
        maxstarnum = maxstar[fitslist0[0][0]]
        threshold_range = [l_threshold, h_threshold, interval]
        starnumlist, coordsfilelist, l_threshold1 = starfind_center3(fitslist0, pixscale[band], satcount[band], threshold_range, minstarnum, maxstarnum)

        #print(f'starnumlist\n{starnumlist}')
        #print(f'coordsflelist\n{coordsfilelist}')

        return starnumlist, coordsfilelist, l_threshold1

    def calc_threshold(fitslist0):
        satcount = {
        'g':param.g_satcount, 'i':param.i_satcount,
        'j':param.j_satcount, 'h':param.h_satcount, 'k':param.k_satcount
        }
        band = fitslist0[0][0]
        stdlist0 = []
        skcount0 = []
        for varr in fitslist0:
            stdlist0.append(bottom.skystat(varr, 'stddev'))
            hdu = fits.open(varr)
            try:
                skcount0.append(float(hdu[0].header['SKYCOUNT']))
            except:
                skcount0.append(0)
        np_stdlist0 = np.array(stdlist0)
        np_skcount0 = np.array(skcount0)
        medstd = np.median(np_stdlist0)
        medskc = np.median(np_skcount0)
        threshold0 = (satcount[band] - medskc)/medstd
        #print(f'{threshold0} = ({satcount[band]} - {medskc})/{medstd}')
        recom_threshold = int(threshold0)

        return recom_threshold

    if optkey:
        for varr in optkey:
            #threshold1 = calc_threshold(fitslist[varr])
            #optstarlist[varr], optcoolist[varr] = iterate_part(fitslist[varr], param, threshold1)
            optstarlist[varr], optcoolist[varr], opt_l_threshold[varr] = iterate_part(fitslist[varr], param, 30, 26, 2)

    if infrakey:
        for varr in infrakey:
            #threshold1 = calc_threshold(fitslist[varr])
            #infstarlist[varr], infcoolist[varr] = iterate_part(fitslist[varr], param, threshold1)
            infstarlist[varr], infcoolist[varr], inf_l_threshold[varr] = iterate_part(fitslist[varr], param, 15, 14)    

    return optstarlist, optcoolist, infstarlist, infcoolist, opt_l_threshold, inf_l_threshold



def check_starnum(optstarlist, optcoolist, infstarlist, infcoolist, opt_l_threshold, inf_l_threshold):
    
    for varr in optstarlist:
        optmed = statistics.median(optstarlist[varr])
        optstd = statistics.stdev(optstarlist[varr])
        optfew = [i for i, num in enumerate(optstarlist[varr]) if optmed - num > 2 * optstd]
        for varr2 in optfew:
            if opt_l_threshold[varr2]==2.0:
                print(f"few stars in {optcoolist[varr][varr2][:-4]}.fits")

    for varr in infstarlist:
        infmed = statistics.median(infstarlist[varr])
        infstd = statistics.stdev(infstarlist[varr])
        inffew = [i for i, num in enumerate(infstarlist[varr]) if infmed - num > 2 * infstd]
        for varr2 in inffew:
            if inf_l_threshold[varr2]==2.0:
                print(f"few stars in {infcoolist[varr][varr2][:-4]}.fits")



def do_xyxymatch(param, optstarlist, optcoolist, infstarlist, infcoolist):

    opt_match = {}
    inf_match = {}
    opt_matchbase = {}
    inf_matchbase = {}
    opt_matchedf = {}
    inf_matchedf = {}

    match_threshold = {
        'g':param.g_threshold, 'i':param.i_threshold,
        'j':param.j_threshold, 'h':param.h_threshold, 'k':param.k_threshold
        }

    if optcoolist:
        optcommon = {s[1:-4] for s in optcoolist[next(iter(optcoolist))]}
        for key in optcoolist:
            optcommon &= {s[1:-4] for s in optcoolist[key]}
        optbase = sorted(optcommon)[0]
    else:
        optcommon = set()

    for varr in optcoolist:
        if optcoolist[varr] and min(optstarlist[varr]) > 3:
            opt_match[varr] = 1
            opt_matchedf[varr] = []
            tempfits = f"{varr}{optbase}.fits"
            hdu = fits.open(tempfits)
            base_rotate = float(hdu[0].header['OFFSETRO']) or 0
            opt_matchbase[varr] = optbase
            for filename in tqdm(optcoolist[varr], desc='{:<}'.format(f'{varr} trimatch')):
                if filename[1:-4] == optbase:
                    continue
                tempfits = re.sub('.coo', '.fits', filename)
                hdu = fits.open(tempfits)
                move_rotate = float(hdu[0].header['OFFSETRO']) or 0
                rotatediff = move_rotate - base_rotate
                outf = re.sub(r'.coo', r'.match', filename)
                referencef = f"{varr}{optbase}.coo"
                outfvarr = triangle_match(filename, referencef, outf, match_threshold[varr])
                if outfvarr == None:
                    continue
                opt_matchedf[varr].append(outfvarr)
        else:
            opt_matchbase[varr] = optbase
            opt_match[varr] = 0


    if infcoolist:
        infcommon = {s[1:-4] for s in infcoolist[next(iter(infcoolist))]}
        for key in infcoolist:
            infcommon &= {s[1:-4] for s in infcoolist[key]}
        infbase = sorted(infcommon)[0]
    else:
        infcommon = set()

    for varr in infcoolist:
        if infcoolist[varr] and min(infstarlist[varr]) > 3:
            inf_match[varr] = 1
            inf_matchedf[varr] = []
            tempfits = f"{varr}{infbase}.fits"
            hdu = fits.open(tempfits)
            base_rotate = float(hdu[0].header['OFFSETRO']) or 0
            inf_matchbase[varr] = infbase
            for filename in tqdm(infcoolist[varr], desc='{:<}'.format(f'{varr} trimatch')):
                if filename[1:-4] == infbase:
                    continue
                tempfits = re.sub('.coo', '.fits', filename)
                hdu = fits.open(tempfits)
                move_rotate = float(hdu[0].header['OFFSETRO']) or 0
                rotatediff = move_rotate - base_rotate
                outf = re.sub(r'.coo', r'.match', filename)
                referencef = f"{varr}{infbase}.coo"
                outfvarr = triangle_match(filename, referencef, outf, match_threshold[varr])
                if outfvarr == None:
                    continue
                inf_matchedf[varr].append(outfvarr)
        else:
            inf_matchbase[varr] = infbase
            inf_match[varr] = 0

    
    if opt_match:
        for varr in opt_matchedf:
            if opt_match[varr] == 1:
                matched_num = match_checker(opt_matchedf[varr])
                if matched_num < 3:
                    opt_match[varr] == 0
                
    if inf_match:
        for varr in inf_matchedf:
            if inf_match[varr] == 1:
                matched_num = match_checker(inf_matchedf[varr])
                if matched_num < 3:
                    inf_match[varr] == 0

    
    return opt_match, opt_matchbase, opt_matchedf, inf_match, inf_matchbase, inf_matchedf


def do_geomap(fitslist, opt_match, opt_matchedf, inf_match, inf_matchedf):

    opt_outlist = {}
    inf_outlist = {}

    for varr in opt_match:
        if opt_match[varr] == 1:
            data = fits.getdata(fitslist[varr][0])
            nxblock = len(data[0])
            nyblock = len(data)
            opt_outlist[varr] = []
            for varrf in opt_matchedf[varr]:
                outf = re.sub(r'.match', r'.geo', varrf)
                bottom.geomap(varrf, nxblock, nyblock, outf, 'rotate')
                opt_outlist[varr].append(outf)
        elif opt_match[varr] == 2:
            data = fits.getdata(fitslist[varr][0])
            nxblock = len(data[0])
            nyblock = len(data)
            opt_outlist[varr] = []
            for varrf in opt_matchedf[varr]:
                outf = re.sub(r'.match', r'.geo', varrf)
                bottom.geomap(varrf, nxblock, nyblock, outf, 'shift')
                opt_outlist[varr].append(outf)
    
    for varr in inf_match:
        if inf_match[varr] == 1:
            data = fits.getdata(fitslist[varr][0])
            nxblock = len(data[0])
            nyblock = len(data)
            inf_outlist[varr] = []
            for varrf in inf_matchedf[varr]:
                outf = re.sub(r'.match', r'.geo', varrf)
                bottom.geomap(varrf, nxblock, nyblock, outf, 'rotate')
                inf_outlist[varr].append(outf)
        elif inf_match[varr] == 2:
            data = fits.getdata(fitslist[varr][0])
            nxblock = len(data[0])
            nyblock = len(data)
            inf_outlist[varr] = []
            for varrf in inf_matchedf[varr]:
                outf = re.sub(r'.match', r'.geo', varrf)
                bottom.geomap(varrf, nxblock, nyblock, outf, 'shift')
                inf_outlist[varr].append(outf)
    
    return opt_outlist, inf_outlist


def do_geotran(fitslist, param, optkey, infrakey, opt_matchb, inf_matchb, opt_geomfile, inf_geomfile):

    def calc_geot(base_band, move_band, coordinate):
        param_name1 = f'{base_band}coo_to_{move_band}'
        param_name2 = f'{move_band}coo_to_{base_band}'
        param_name3 = f'theta_{base_band}{move_band}'
        param_name4 = f'theta_{move_band}{base_band}'
        baseb_coo = getattr(param, param_name1)
        moveb_coo = getattr(param, param_name2)
        baseb_coox = float(baseb_coo[0])
        baseb_cooy = float(baseb_coo[1])
        moveb_coox = float(moveb_coo[0])
        moveb_cooy = float(moveb_coo[1])
        
        try:
            temp_theta = getattr(param, param_name3)
            theta_degree = float(temp_theta)
        except:
            temp_theta = getattr(param, param_name4)
            theta_degree = -float(temp_theta)
            
        theta_rad = np.radians(theta_degree)
        x    = coordinate[0]
        y    = coordinate[1]
        x_in = baseb_coox
        y_in = baseb_cooy
        moved_coox = (x - x_in)*np.cos(theta_rad)-(y - y_in)*np.sin(theta_rad)+moveb_coox
        moved_cooy = (x - x_in)*np.sin(theta_rad)+(y - y_in)*np.cos(theta_rad)+moveb_cooy

        return moved_coox, moved_cooy

    opt_geomdict = {}
    
    for varr in opt_geomfile:
        tempfits = f'{varr}{opt_matchb[varr]}.fits'
        hdu = fits.open(tempfits)
        base_rotate = float(hdu[0].header['OFFSETRO']) or 0
        opt_geomdict[varr] = geotparam(param, opt_geomfile[varr], base_rotate)
    inf_geomdict = {}
    for varr in inf_geomfile:
        tempfits = f'{varr}{inf_matchb[varr]}.fits'
        hdu = fits.open(tempfits)
        base_rotate = float(hdu[0].header['OFFSETRO']) or 0
        inf_geomdict[varr] = geotparam(param, inf_geomfile[varr], base_rotate)
    

    opt_iddict = {}
    for varr in opt_geomdict:
        for index, varr2 in enumerate(opt_geomdict[varr]):
            if varr2['fitsid'] not in opt_iddict:
                opt_iddict[varr2['fitsid']] = {}
            opt_iddict[varr2['fitsid']][varr] = index

    inf_iddict = {}
    for varr in inf_geomdict:
        #print(f'え？{inf_geomdict}')
        for index, varr2 in enumerate(inf_geomdict[varr]):
            #print(f'まじ？{varr2}')
            if varr2['fitsid'] not in inf_iddict:
                inf_iddict[varr2['fitsid']] = {}
            inf_iddict[varr2['fitsid']][varr] = index


    geotran_base = {
        'g':param.tran_g, 'i':param.tran_i,
        'j':param.tran_j, 'h':param.tran_h, 'k':param.tran_k
    }

    bottom.geotran_param()

    not_exec = []
    basefits = {}
    for varr in optkey:
        
        for fitsname in tqdm(fitslist[varr], desc=f'try {varr} band geotran '):
            fitsid = fitsname[1:-5]
            
            if fitsid == opt_matchb[geotran_base[varr]]:
                outfile = re.sub('.fits', f'_geo{varr}.fits', fitsname)
                bottom.geotran(fitsname, outfile, 1, 1, 1, 1, 0, 0)
                basefits[varr] = outfile
                continue

            if fitsid not in opt_iddict:
                not_exec.append(fitsname)
                continue
            
            if geotran_base[varr] in opt_iddict[fitsid]:
                base_band = geotran_base[varr]
                index = opt_iddict[fitsid][base_band]
            else:
                not_exec.append(fitsname)
                continue

            outfile = re.sub('.fits', f'_geo{base_band}.fits', fitsname)
            xmean = opt_geomdict[base_band][index]['xmean']
            ymean = opt_geomdict[base_band][index]['ymean']
            xrefmean = opt_geomdict[base_band][index]['xrefmean']
            yrefmean = opt_geomdict[base_band][index]['yrefmean']
            xrotation = opt_geomdict[base_band][index]['xrotation']
            yrotation = opt_geomdict[base_band][index]['yrotation']

            if varr != base_band:
                xmean, ymean = calc_geot(base_band, varr, (xmean, ymean))
                xrefmean, yrefmean = calc_geot(base_band, varr, (xrefmean, yrefmean))

            bottom.geotran(fitsname, outfile, xmean, ymean, xrefmean, yrefmean, xrotation, yrotation)

    #print(f'inf_iddict\n{inf_iddict}')
    #print(f'inf_geomdict\n{inf_geomdict}')
    for varr in infrakey:
        
        for fitsname in tqdm(fitslist[varr], desc=f'try {varr} band geotran '):
            fitsid = fitsname[1:-5]

            if fitsid == inf_matchb[geotran_base[varr]]:
                outfile = re.sub('.fits', f'_geo{varr}.fits', fitsname)
                bottom.geotran(fitsname, outfile, 1, 1, 1, 1, 0, 0)
                basefits[varr] = outfile
                continue

            if fitsid not in inf_iddict:
                not_exec.append(fitsname)
                continue

            if geotran_base[varr] in inf_iddict[fitsid]:
                base_band = geotran_base[varr]
                index = inf_iddict[fitsid][base_band]
            elif inf_iddict[fitsid]:
                #print(fitsname)
                base_band = list(inf_iddict[fitsid].keys())[0]
                index = inf_iddict[fitsid][base_band]
            else:
                not_exec.append(fitsname)
                continue

            outfile = re.sub('.fits', f'_geo{base_band}.fits', fitsname)
            xmean = inf_geomdict[base_band][index]['xmean']
            ymean = inf_geomdict[base_band][index]['ymean']
            xrefmean = inf_geomdict[base_band][index]['xrefmean']
            yrefmean = inf_geomdict[base_band][index]['yrefmean']
            xrotation = inf_geomdict[base_band][index]['xrotation']
            yrotation = inf_geomdict[base_band][index]['yrotation']

            if varr != base_band:
                xmean, ymean = calc_geot(base_band, varr, (xmean, ymean))
                xrefmean, yrefmean = calc_geot(base_band, varr, (xrefmean, yrefmean))
                
            bottom.geotran(fitsname, outfile, xmean, ymean, xrefmean, yrefmean, xrotation, yrotation)

    
    for band in basefits:
        print(f'{band} base is {basefits[band]}')
    
    not_exec.sort()
    for varr in not_exec:
        print(f'{varr} was not moved.')
    


def main(fitslist, param):
     
    subprocess.run('rm *.match' ,shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm *geo*' ,shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm *.coo' ,shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm *_xysh*' ,shell=True, stderr=subprocess.DEVNULL)

    #base の認識

    keyslist = list(fitslist.keys())
    optkey = list(set(keyslist) & set(['g', 'i']))
    infrakey = list(set(keyslist) & set(['j', 'h', 'k']))
    
    result_varr = do_starfind(fitslist, param, optkey, infrakey)

    optstarlist = result_varr[0]
    optcoolist  = result_varr[1]
    infstarlist = result_varr[2]
    infcoolist  = result_varr[3]

    opt_l_threshold = result_varr[4]
    inf_l_threshold = result_varr[5]

    #星が見つからないとき、この辺で終了した方がいいかもな。
    if not optcoolist and not infcoolist:
        print('if not optcoolist and not infracoolist:')
        sys.exit()

    check_starnum(optstarlist, optcoolist, infstarlist, infcoolist, opt_l_threshold, inf_l_threshold)

    
    result_varr = do_xyxymatch(param, optstarlist, optcoolist, infstarlist, infcoolist)

    opt_match    = result_varr[0]
    opt_matchb   = result_varr[1]
    opt_matchedf = result_varr[2]
    inf_match    = result_varr[3]
    inf_matchb   = result_varr[4]
    inf_matchedf = result_varr[5]

    if all(value == 0 for value in opt_match.values()) and all(value == 0 for value in inf_match.values()):
        print(f'all matches failed')
        sys.exit()

    #何をreturn するかはその後ができてから。
    result_varr = do_geomap(fitslist, opt_match, opt_matchedf, inf_match, inf_matchedf)

    opt_geomfile = result_varr[0]
    inf_geomfile = result_varr[1]
    
    do_geotran(fitslist, param, optkey, infrakey, opt_matchb, inf_matchb, opt_geomfile, inf_geomfile)

    