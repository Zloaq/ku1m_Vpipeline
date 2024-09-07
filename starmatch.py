#!/opt/anaconda3/envs/p11/bin/python3

import os
import re
import sys
import math
import subprocess
import numpy as np
import pandas as pd
from scipy import ndimage
from tqdm import tqdm
from astropy.io import fits
from pyraf import iraf

import bottom

import matplotlib.pyplot as plt

def starfind_center3(fitslist, param, searchrange=[3.0, 5.0, 0.2], minstarnum=0, maxstarnum=100):
    
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
        tempkey = f'{band}_satcount'
        #div には対応してない
        skycount = float(header['SKYCOUNT'] or 0)
        saturation_value = float(getattr(param, tempkey)) - skycount
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

    pixscale = {
        'g':param.pixscale_g, 'i':param.pixscale_i,
        'j':param.pixscale_j, 'h':param.pixscale_h, 'k':param.pixscale_k
    }
    coordsfilelist = []
    starnumlist = []
    band = fitslist[0][0]
    #print(f'{band} band threshold range = {searchrange}')

    #for index, filename in enumerate(tqdm(fitslist, desc=f'{fitslist[0][0]} band starfind')):
    iterate = 0
    searchrange0 = searchrange
    with tqdm(total=len(fitslist), desc=f'{fitslist[0][0]} band starfind') as pbar:
        while iterate <= len(fitslist) - 1:
            filename = fitslist[iterate]
            data = fits.getdata(filename)
            header = fits.getheader(filename)
            offra_pix = int(float(header['OFFSETRA'])/pixscale[band])
            offde_pix = int(float(header['OFFSETDE'])/pixscale[band])
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
            pbar.update(1)

    return starnumlist, coordsfilelist, iterate


def lnum_match(inputf, referancef, outputf, rotatediff=0): # 変な画像があることは考慮してない。

    def read_coofile(infile):
        with open(infile, 'r') as file:
            flines = file.readlines()
        lines = [line.strip().split() for line in flines if not line.startswith('#')]
        coords = [[float(line[0]), float(line[1])] for line in lines] # coox, cooy,
        return coords

    
    def calculate_segments(coords):
        coords = np.array(coords)
        points = coords.shape[0]

        if points <= 1:
            return np.array([])

        p1_idx, p2_idx = np.triu_indices(points, k=1)

        x1, y1 = coords[p1_idx, 0], coords[p1_idx, 1]
        x2, y2 = coords[p2_idx, 0], coords[p2_idx, 1]

        lengths = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        angle_rads = np.arctan2(y2 - y1, x2 - x1)
        angles_deg = np.degrees(angle_rads)

        
        mask = (y1 < y2) | ((y1 == y2) & (x1 < x2))

        segments = np.zeros((len(p1_idx), 3, 2))
        segments[mask, 0] = coords[p1_idx[mask]]
        segments[mask, 1] = coords[p2_idx[mask]]
        segments[~mask, 0] = coords[p2_idx[~mask]]
        segments[~mask, 1] = coords[p1_idx[~mask]]
        
        angle1 = np.arctan2(segments[:,1,1]-segments[:,0,1], segments[:,1,0]-segments[:,0,0])
        angle2 = np.arctan2(segments[:,1,0]-segments[:,0,0], segments[:,1,1]-segments[:,0,1])
        angle_mask = (np.pi/4 < angle1) & (angle1 < 3*np.pi/4)

        segments[angle_mask, 2, 1] = np.degrees(angle1[angle_mask])
        segments[~angle_mask, 2, 1] = 90 - np.degrees(angle2[~angle_mask])

        segments[:, 2, 0] = lengths
        segments[:, 2, 1] = angles_deg

        return segments
    

    def match_segments(inlist, reflist, length_tolerance=5e-1, angle_tolerance=1e-1):

        inplist = np.array(inlist)
        reflist = np.array(reflist)
        in_lengths = inplist[:, 2, 0]
        ref_lengths = reflist[:, 2, 0]
        in_angles = inplist[:, 2, 1]
        ref_angles = reflist[:, 2, 1]

        length_diffs = np.abs(ref_lengths[:, np.newaxis] - in_lengths[np.newaxis, :])
        
        length_mask = length_diffs < length_tolerance

        indices1 = np.argwhere(length_mask)


        #長さが同じものの位置つまり長さが同じ組み合わせ

        reflist1 = reflist[indices1[:, 0], :]
        inplist1 = inplist[indices1[:, 1], :]

        refvect = reflist1[:, 1] - reflist1[:, 0]
        inpvect = inplist1[:, 1] - inplist1[:, 0]
        reflengths = reflist1[:, 2, 0][:, np.newaxis]
        inplengths = inplist1[:, 2, 0][:, np.newaxis]
        refvect_norm = refvect / reflengths
        inpvect_norm = inpvect / inplengths

        refcos = np.dot(refvect_norm, refvect_norm.T)
        inpcos = np.dot(inpvect_norm, inpvect_norm.T)
        mask0 = np.triu(np.ones(refcos.shape, dtype=bool), k=1)
        indices2 = np.argwhere(mask0)
        #indices2 同じ長さの組の内積

        refcos1 = refcos[mask0]
        inpcos1 = inpcos[mask0]

        #cos_mask = np.abs(refcos1) - np.abs(inpcos1) < 1e-1

        cos_mask0 = np.abs(refcos1) - np.abs(inpcos1) < 1e-1
        cos_mask1 = np.abs(refcos1)!=1
        cos_mask = cos_mask0 & cos_mask1

        refcos2 = refcos1[cos_mask]
        inpcos2 = inpcos1[cos_mask]
        indices3 = indices2[cos_mask]
        #indices3 内積計算したやつの生き残り
        refseg0 = reflist1[indices3[:, 0], :]
        refseg1 = reflist1[indices3[:, 1], :]
        inpseg0 = inplist1[indices3[:, 0], :]
        inpseg1 = inplist1[indices3[:, 1], :]

        refdistance1 = np.linalg.norm(
            refseg0[:, [0, 0, 1, 1], :] - refseg1[:, [0, 1, 0, 1], :], axis=2
            )
        inpdistance1 = np.linalg.norm(
            inpseg0[:, [0, 0, 1, 1], :] - inpseg1[:, [0, 1, 0, 1], :], axis=2
            )
        
        refdistance2 = np.sort(refdistance1, axis=1)
        inpdistance2 = np.sort(inpdistance1, axis=1)
        ref_indices2 = np.argsort(refdistance1, axis=1)
        inp_indices2 = np.argsort(inpdistance1, axis=1)

        mask1 = np.abs(refdistance2 - inpdistance2) < length_tolerance
        mask2 = mask1.all(axis=1)
        indices4 = indices3[mask2]
        
        ref_indices3 = np.empty((2*indices4.shape[0], 4))
        inp_indices3 = np.empty((2*indices4.shape[0], 4))
        ref_indices3[0::2] = ref_indices2[indices4[:, 0]]
        ref_indices3[1::2] = ref_indices2[indices4[:, 1]]
        inp_indices3[0::2] = inp_indices2[indices4[:, 0]]
        inp_indices3[1::2] = inp_indices2[indices4[:, 1]]
        #これはrefcos = np.dot(refvect_norm, refvect_norm.T)inpcos = np.dot(inpvect_norm, inpvect_norm.T))でした。
        
        reflist2 = np.empty((2*indices4.shape[0], 3, 2))
        inplist2 = np.empty((2*indices4.shape[0], 3, 2))
        reflist2[0::2] = reflist1[indices4[:, 0]]
        reflist2[1::2] = reflist1[indices4[:, 1]]
        inplist2[0::2] = inplist1[indices4[:, 0]]
        inplist2[1::2] = inplist1[indices4[:, 1]]

        first_element_diff = ref_indices3[:, 0] - inp_indices3[:, 0]
        mask3 = first_element_diff == 0


        reflist3 = reflist2[:, :2, :]
        inplist3 = np.empty_like(reflist3)

        inplist3[mask3, :, :] = inplist2[mask3][:, :2, :]
        inplist3[~mask3, :, :] = inplist2[~mask3][:, [1, 0], :]

        reflist4 = reflist3[:, 0, :]
        inplist4 = inplist3[:, 0, :]
        reflist5 = reflist3[:, 1, :]
        inplist5 = inplist3[:, 1, :]

        matched_pair1 = np.array([reflist4, inplist4])
        matched_pair2 = np.array([reflist5, inplist5])

        matched_pair1 = matched_pair1.transpose(1, 0, 2)
        matched_pair2 = matched_pair2.transpose(1, 0, 2)

        matched_pair = np.empty((matched_pair1.shape[0] + matched_pair2.shape[0], 2, 2), dtype=matched_pair1.dtype)
        matched_pair[0::2] = matched_pair1
        matched_pair[1::2] = matched_pair2


        return matched_pair

        
    def one_point_match(incoords, refcoords):#一番明るさの近いやつでいいや 距離？じゃね？
        #print('one point match')
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
        return None#えーん
        matched_coo = one_point_match(inpcoords, refcoords)
    else:
        matched_coo = match_segments(inpseg, refseg)
    # write .match

    if matched_coo.size == 0:
        print()
        return None

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


def do_lnummatch(optcoolist, infcoolist):
    for varr in optcoolist:
        for filename in optcoolist[varr][1:]:
            outf = re.sub(r'.coo', r'.match', filename)
            lnum_match(filename, optcoolist[varr][0], outf)


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
    """
    def iterate_part(fitslist0, param, h_threshold=10, l_threshold=9, interval=0.2):
        starnumlist = [0]
        minstarnum = 3
        maxstarnum = 40
        renum = [0, 0]
        l_threshold_org = l_threshold
        while (min(starnumlist) < minstarnum or max(starnumlist) > maxstarnum) and l_threshold > 2.0:
            threshold_range = [l_threshold, h_threshold, interval]
            starnumlist, coordsfilelist, index0 = starfind_center3(fitslist0, param, threshold_range, minstarnum, maxstarnum)
            #print(starnumlist)###
            fitslist0 = fitslist0[index0:] + fitslist0[:index0]

            if min(starnumlist) < minstarnum:
                l_threshold -= 1
                h_threshold -= 1
                renum[0] = 1
                if renum == [1, 1]:
                    print(f'{fitslist0[0]} is low')
                    del fitslist0[0]
                    print(f'reject {fitslist0[0]}')       

            elif max(starnumlist) > maxstarnum:
                l_threshold += 1
                h_threshold += 1
                renum[1] = 1
                if renum == [1, 1]:
                    print(f'{fitslist0[0]} is high')

        return starnumlist, coordsfilelist
    """

    def iterate_part(fitslist0, param, h_threshold=10, l_threshold=9, interval=0.2):
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
        starnumlist, coordsfilelist, index0 = starfind_center3(fitslist0, param, threshold_range, minstarnum, maxstarnum)

        #print(f'starnumlist\n{starnumlist}')
        #print(f'coordsflelist\n{coordsfilelist}')

        return starnumlist, coordsfilelist

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
            optstarlist[varr], optcoolist[varr] = iterate_part(fitslist[varr], param, 30, 26, 2)

    if infrakey:
        for varr in infrakey:
            #threshold1 = calc_threshold(fitslist[varr])
            #infstarlist[varr], infcoolist[varr] = iterate_part(fitslist[varr], param, threshold1)
            infstarlist[varr], infcoolist[varr] = iterate_part(fitslist[varr], param, 15, 14)    

    return optstarlist, optcoolist, infstarlist, infcoolist


def do_xyxymatch(optstarlist, optcoolist, infstarlist, infcoolist):

    opt_match = {}
    inf_match = {}
    opt_matchbase = {}
    inf_matchbase = {}
    opt_matchedf = {}
    inf_matchedf = {}

    if optcoolist:
        optcommon = set(s[1:-4] for s in optcoolist[next(iter(optcoolist))])
        optbase   = sorted(optcommon)[0]
    else:
        optcommon = set()

    

    for varr in optcoolist:
        if optcoolist[varr] and min(optstarlist[varr]) > 3:
            opt_match[varr] = 1
            opt_matchedf[varr] = []
            tempfits = re.sub('.coo', '.fits', optcoolist[varr][0])
            hdu = fits.open(tempfits)
            base_rotate = float(hdu[0].header['OFFSETRO']) or 0
            opt_matchbase[varr] = optbase
            for filename in tqdm(optcoolist[varr], desc='{:<}'.format(f'{varr} lnummatch')):
                if filename[1:-4] == optbase:
                    continue
                tempfits = re.sub('.coo', '.fits', filename)
                hdu = fits.open(tempfits)
                move_rotate = float(hdu[0].header['OFFSETRO']) or 0
                rotatediff = move_rotate - base_rotate
                outf = re.sub(r'.coo', r'.match', filename)
                outfvarr = lnum_match(filename, optcoolist[varr][0], outf, rotatediff)
                if outfvarr == None:
                    continue
                opt_matchedf[varr].append(outfvarr)
        else:
            opt_matchbase[varr] = optcommon
            opt_match[varr] = 0

    if infcoolist:
        infcommon = set(s[1:-4] for s in infcoolist[next(iter(infcoolist))])
        infbase   = sorted(infcommon)[0]
    else:
        infcommon = set()

    
    #print(f'infcommon = {infcommon}')

    for varr in infcoolist:
        if infcoolist[varr] and min(infstarlist[varr]) > 3:
            inf_match[varr] = 1
            inf_matchedf[varr] = []
            tempfits = re.sub('.coo', '.fits', infcoolist[varr][0])
            hdu = fits.open(tempfits)
            base_rotate = float(hdu[0].header['OFFSETRO']) or 0
            inf_matchbase[varr] = infbase
            for filename in tqdm(infcoolist[varr], desc='{:<}'.format(f'{varr} lnummatch')):
                if filename[1:-4] == infbase:
                    continue
                tempfits = re.sub('.coo', '.fits', filename)
                hdu = fits.open(tempfits)
                move_rotate = float(hdu[0].header['OFFSETRO']) or 0
                rotatediff = move_rotate - base_rotate
                outf = re.sub(r'.coo', r'.match', filename)
                outfvarr = lnum_match(filename, infcoolist[varr][0], outf, rotatediff)
                if outfvarr == None:
                    continue
                inf_matchedf[varr].append(outfvarr)
        else:
            inf_matchbase[varr] = infcommon
            inf_match[varr] = 0

    #この辺で match checker を入れて、coo ファイルにちゃんと書き込まれているか確認
    #星検出してるのに、match できていないものを検出。
    #coolist に入っているもののうち、matchが行われていないものと、matchできていないもの。
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
    for varr in optkey:
        #ほんとに 1からでいいんか？？？？
        for fitsname in tqdm(fitslist[varr][1:], desc=f'{varr} band geotran '):
            fitsid = fitsname[1:-5]
            
            if fitsid == opt_matchb[geotran_base[varr]]:
                outfile = re.sub('.fits', f'_geo{varr}.fits', fitsname)
                bottom.geotran(fitsname, outfile, 1, 1, 1, 1, 0, 0)
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
        #ほんとに 1からでいいんか？？？？
        for fitsname in tqdm(fitslist[varr][1:], desc=f'{varr} band geotran '):
            fitsid = fitsname[1:-5]

            if fitsid == inf_matchb[geotran_base[varr]]:
                outfile = re.sub('.fits', f'_geo{varr}.fits', fitsname)
                bottom.geotran(fitsname, outfile, 1, 1, 1, 1, 0, 0)
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

    not_exec.sort()
    for varr in not_exec:
        print(f'{varr} was not moved.')
    


    

def main(fitslist, param):
     
    subprocess.run('rm *.match' ,shell=True)
    subprocess.run('rm *geo*' ,shell=True)
    subprocess.run('rm *.coo' ,shell=True)
    subprocess.run('rm *_xysh*' ,shell=True)

    #base の認識

    keyslist = list(fitslist.keys())
    optkey = list(set(keyslist) & set(['g', 'i']))
    infrakey = list(set(keyslist) & set(['j', 'h', 'k']))
    
    result_varr = do_starfind(fitslist, param, optkey, infrakey)

    optstarlist = result_varr[0]
    optcoolist  = result_varr[1]
    infstarlist = result_varr[2]
    infcoolist  = result_varr[3]

    #星が見つからないとき、この辺で終了した方がいいかもな。
    if not optcoolist and not infcoolist:
        print('if not optcoolist and not infracoolist:')
        sys.exit()

    #print(f'starnum \n{infstarlist}')
    
    result_varr = do_xyxymatch(optstarlist, optcoolist, infstarlist, infcoolist)

    opt_match    = result_varr[0]
    opt_matchb   = result_varr[1]
    opt_matchedf = result_varr[2]
    inf_match    = result_varr[3]
    inf_matchb   = result_varr[4]
    inf_matchedf = result_varr[5]

    #何をreturn するかはその後ができてから。
    result_varr = do_geomap(fitslist, opt_match, opt_matchedf, inf_match, inf_matchedf)

    opt_geomfile = result_varr[0]
    inf_geomfile = result_varr[1]
    
    do_geotran(fitslist, param, optkey, infrakey, opt_matchb, inf_matchb, opt_geomfile, inf_geomfile)

    