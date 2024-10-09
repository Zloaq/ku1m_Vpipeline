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

def starfind_center3(fitslist, pixscale, satcount, searchrange=[3.0, 5.0, 0.2], minstarnum=4, maxstarnum=100, minthreshold=2.0):
    
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
        return np.where(data > threshold * rms, data, 0)
    
    def binarize_data(data, threshold, rms):
        return data > threshold * rms
    
    def filter_saturate(labeled_image, filtered_data, object_slices, header):
        #div には対応してない
        skycount = float(header['SKYCOUNT'] or 0)
        saturation_value = float(satcount) - skycount
        saturated_mask = filtered_data >= saturation_value
        saturated_labels = np.unique(labeled_image[saturated_mask])
        
        templabels = np.arange(1, len(object_slices) + 1)
        mask = ~np.isin(templabels, saturated_labels)
        filtered_labels = templabels[mask]
        filtered_objects = np.array(object_slices)[mask]
        filtered_objects = [tuple(inner_list) for inner_list in filtered_objects.tolist()]
        #print(f'filtered_labels\n{filtered_labels.tolist()}')
        #print(f'filtered_objects\n{tuple(filtered_objects.tolist())}')
        #sys.exit()
        
        return filtered_labels.tolist(), filtered_objects

    def detect_round_clusters(filtered_labels, filtered_objects, labeled_image, square=5, fillrate=0.25):
        squareness_values = np.array([squareness(region_slice) for region_slice in filtered_objects])
        filling_rates = np.array([filling_rate(filtered_labels[i], region_slice, labeled_image) 
                                for i, region_slice in enumerate(filtered_objects)])

        mask = (squareness_values < square) & (filling_rates > fillrate)
        
        round_clusters = np.array(filtered_labels)[mask]
        slice_list = np.array([[region_slice[0], region_slice[1]] for region_slice in filtered_objects])[mask]

        return round_clusters, slice_list
    
    def clustar_centroid(data, slices, padding=2):
        max_y, max_x = data.shape

        # スライスのリストを numpy 配列に変換
        slices_array = np.array([[sl[0].start, sl[0].stop, sl[1].start, sl[1].stop] for sl in slices])
        #print(f'aaaa{slices_array}')
        # パディング適用後のスライス範囲を計算
        y_starts = np.clip(slices_array[:, 0] - padding, 0, max_y)
        y_stops = np.clip(slices_array[:, 1] + padding, 0, max_y)
        x_starts = np.clip(slices_array[:, 2] - padding, 0, max_x)
        x_stops = np.clip(slices_array[:, 3] + padding, 0, max_x)

        centroids = []

        for y_start, y_stop, x_start, x_stop in zip(y_starts, y_stops, x_starts, x_stops):
            data_slice = data[y_start:y_stop, x_start:x_stop]
            
            total = data_slice.sum()
            if total == 0:
                continue

            y_indices, x_indices = np.indices(data_slice.shape)
            y_centroid_local = np.sum(y_indices * data_slice) / total
            x_centroid_local = np.sum(x_indices * data_slice) / total
            y_centroid_global = y_centroid_local + y_start
            x_centroid_global = x_centroid_local + x_start
            centroids.append((y_centroid_global, x_centroid_global))

        return centroids

    """
    def clustar_centroid(data, slices, padding=2):
        centroids = []
        for sl in slices:
            y_slice, x_slice = sl
            y_start = max(y_slice.start - padding, 0)
            y_stop = min(y_slice.stop + padding, data.shape[0])
            x_start = max(x_slice.start - padding, 0)
            x_stop = min(x_slice.stop + padding, data.shape[1])
            data_slice = data[y_start:y_stop, x_start:x_stop]
            y_indices, x_indices = np.indices(data_slice.shape)
            total = np.sum(data_slice)
            if total == 0:
                continue
            y_centroid_local = np.sum(y_indices * data_slice) / total
            x_centroid_local = np.sum(x_indices * data_slice) / total
            y_centroid_global = y_centroid_local + y_slice.start
            x_centroid_global = x_centroid_local + x_slice.start
            centroids.append((y_centroid_global, x_centroid_global))
        
        return centroids
    """

    def chose_unique_coords(center_list, threshold=3):
        # numpy配列に変換
        center_array = np.array(center_list)
        unique_centers = []

        while len(center_array) > 0:
            ref_point = center_array[0]
            unique_centers.append(ref_point)

            distances = np.sqrt(np.sum((center_array - ref_point) ** 2, axis=1))

            center_array = center_array[distances > threshold]

        return unique_centers
    
    def write_to_txt(centers, filename):
        with open(filename, 'w') as f1:
            for center in centers:
                f1.write(f'{center[1]}  {center[0]}\n')
    
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
            #print(f'filenamr {filename}')
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

            loopnum = [0, 0]
            while searchrange0[0] >= minthreshold:
                center_list = []
                #print()
                for threshold in np.arange(searchrange0[0], searchrange0[1], searchrange0[2]):
                    binarized_data = binarize_data(data, threshold, rms)
                    labeled_image, _ = ndimage.label(binarized_data)
                    object_slices = ndimage.find_objects(labeled_image)
                    filtered_data = filter_data(data, threshold, rms)
                    filtered_labels, filtered_objects = filter_saturate(labeled_image, filtered_data, object_slices, header)

                    _, slice_list = detect_round_clusters(filtered_labels, filtered_objects, labeled_image)

                    if len(slice_list) == 0:
                        continue

                    centers = clustar_centroid(data, slice_list)
                    center_list.extend(centers)

                #print(len(center_list))
                

                unique_center_list = chose_unique_coords(center_list)
                starnum = len(unique_center_list)
                #print(f'{filename}, {searchrange0}')

                if loopnum[0] > 0 and loopnum[1] > 0:
                    if searchrange0[0] < minthreshold:
                        searchrange0[0] -= 0.5
                        searchrange0[1] -= 0.5
                    break

                elif (starnum < minstarnum) and (searchrange0[0] > minthreshold):
                    #print(f'retry starfind')
                    searchrange0[0] -= 0.5
                    searchrange0[1] -= 0.5
                    loopnum[0] += 1
                    continue
                elif starnum > maxstarnum:
                    searchrange0[0] += 0.5
                    searchrange0[1] += 0.5
                    loopnum[1] += 1
                    continue
                else:
                    break
            
            file = f'{filename[:-5]}.coo'
            write_to_txt(unique_center_list, file)
            coordsfilelist.append(file)
            starnumlist.append(starnum)
            iterate += 1
            threshold_lside.append(searchrange0[0])
            pbar.update(1)

    return starnumlist, coordsfilelist, threshold_lside


def triangle_match(inputf, referencef, outputf, match_threshold=0.05, shift_threshold=5, rotate_threshold=0.5):
    def read_coofile(infile):
        with open(infile, 'r') as file:
            flines = file.readlines()
        lines = [line.strip().split() for line in flines if not line.startswith('#')]
        coords = np.array([[float(line[0]), float(line[1])] for line in lines])  # coox, cooy,
        return coords

    def compute_triangle_descriptors(coords, eps=4):
        triangles = list(combinations(coords, 3))
        descriptors = []

        for tri in triangles:
            #print(0)
            tri = np.array(tri)
            # 辺の長さを計算
            sides = [np.linalg.norm(tri[i] - tri[j]) for i, j in combinations(range(3), 2)]
            #print(1)
            sides = np.sort(sides)
            if any(side < eps for side in sides):
                continue
            a, b, c = sides
            angles = [
            np.arccos(np.clip((b**2 + c**2 - a**2) / (2 * b * c), -1, 1)),
            np.arccos(np.clip((a**2 + c**2 - b**2) / (2 * a * c), -1, 1)),
            np.arccos(np.clip((a**2 + b**2 - c**2) / (2 * a * b), -1, 1))
            ]
            #print(3)
            angles = np.sort(angles)
            #print(4)
            descriptor = np.concatenate([sides, angles])
            descriptors.append((descriptor, tri))
        #print(f'stop')
        return descriptors
    
    def estimate_affine_matrices(src_triangles, dst_triangles):
        affine_matrices = []
        for src, dst in zip(src_triangles, dst_triangles):
            A = []
            b = []
            for (x, y), (x_p, y_p) in zip(src, dst):
                A.append([x, y, 1, 0, 0, 0])
                A.append([0, 0, 0, x, y, 1])
                b.append(x_p)
                b.append(y_p)
            A = np.array(A)
            b = np.array(b)
            h, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            affine_matrix = np.array([[h[0], h[1], h[2]],
                                      [h[3], h[4], h[5]]])
            affine_matrices.append(affine_matrix)
        return np.array(affine_matrices)
        

    def evaluate_affine_matrix(matrix, src_points, dst_points, threshold):
        transformed_points = (matrix[:, :2] @ src_points.T).T + matrix[:, 2]
        distances = np.linalg.norm(transformed_points - dst_points, axis=1)
        return np.sum(distances < threshold)

    def calculate_mode(data, binsize):
        bin_edges = np.arange(np.min(data), np.max(data) + 2*binsize, binsize)
        counts, bin_edges = np.histogram(data, bins=bin_edges)
        #print(f'data {data}')
        #print(f'bin {bin_edges}')
        print(f'binsize {binsize}')
        print(f'counts {counts}')
        max_bin_index = np.argmax(counts)
        mode_value = (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]) / 2
        return mode_value
    
    def match_triangle():
        coords_input = read_coofile(inputf)
        coords_ref = read_coofile(referencef)

        descriptors_input = compute_triangle_descriptors(coords_input)
        descriptors_ref = compute_triangle_descriptors(coords_ref)

        descs_input = np.array([desc[0] for desc in descriptors_input])
        descs_ref = np.array([desc[0] for desc in descriptors_ref])


        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(descs_ref)
        distances, indices = nbrs.kneighbors(descs_input)

        good_matches = distances[:, 0] < match_threshold

        matched_triangles_input = np.array([descriptors_input[i][1] for i in range(len(descriptors_input)) if good_matches[i]])
        matched_triangles_ref = np.array([descriptors_ref[indices[i][0]][1] for i in range(len(descriptors_input)) if good_matches[i]])

        if len(matched_triangles_input) == 0:
            #print(inputf)
            return None

        #print(f'len inp {len(matched_triangles_input)} ref {len(matched_triangles_ref)}')
        # アフィン変換行列を一括して計算
        affine_matrices = estimate_affine_matrices(matched_triangles_input, matched_triangles_ref)

        #print(f'affin {affine_matrices}')

        # アフィン変換行列の評価と選択 (成分(1, 3), (2, 3)および回転成分の最頻値を使用して外れ値を除去)
        t_x = affine_matrices[:, 0, 2]
        t_y = affine_matrices[:, 1, 2]
        rotation_cos = affine_matrices[:, 0, 0]
        rotation_sin = affine_matrices[:, 1, 0]

        # 最頻値を計算
        #base_t_x = calculate_mode(t_x, binsize=5)
        #base_t_y = calculate_mode(t_y, binsize=5)
        #base_rotation_cos = calculate_mode(rotation_cos, binsize=0.2)
        #base_rotation_sin = calculate_mode(rotation_sin, binsize=0.2)

        # 中央値でもいいかもしれない

        base_t_x = np.median(t_x)
        base_t_y = np.median(t_y)
        base_rotation_cos = np.median(rotation_cos)
        base_rotation_sin = np.median(rotation_sin)
        


        # 最頻値からの差が大きいものを外れ値として除去
        valid_indices = (
            (np.abs(t_x - base_t_x) < shift_threshold) &
            (np.abs(t_y - base_t_y) < shift_threshold) &
            (np.abs(rotation_cos - base_rotation_cos) < rotate_threshold) &
            (np.abs(rotation_sin - base_rotation_sin) < rotate_threshold)
        )
        filtered_affine_matrices = affine_matrices[valid_indices]
        filtered_triangles_input = matched_triangles_input[valid_indices]
        filtered_triangles_ref = matched_triangles_ref[valid_indices]

        if len(filtered_triangles_input) == 0:
            return None

        src_points = []
        dst_points = []
        for tri_input, tri_ref in zip(filtered_triangles_input, filtered_triangles_ref):
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

    if matched_coo is None:
        #print(inputf)
        return None

    if matched_coo.size == 0:
        return None
    
    #print(f'coo {matched_coo}')

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

    def iterate_part(fitslist0, param, l_threshold=9, h_threshold=10, interval=0.2):
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
        minstarnum = minstar[fitslist0[0][0]]
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
            optstarlist[varr], optcoolist[varr], opt_l_threshold[varr] = iterate_part(fitslist[varr], param, 26, 30, 2)

    if infrakey:
        for varr in infrakey:
            #threshold1 = calc_threshold(fitslist[varr])
            #infstarlist[varr], infcoolist[varr] = iterate_part(fitslist[varr], param, threshold1)
            infstarlist[varr], infcoolist[varr], inf_l_threshold[varr] = iterate_part(fitslist[varr], param, 14, 15, 0.2)    

    return optstarlist, optcoolist, infstarlist, infcoolist, opt_l_threshold, inf_l_threshold



def check_starnum(optstarlist, optcoolist, infstarlist, infcoolist, opt_l_threshold, inf_l_threshold, minthreshold=2.0):
    
    for varr in optstarlist:
        optmed = statistics.median(optstarlist[varr])
        optstd = statistics.stdev(optstarlist[varr])
        optfew = [i for i, num in enumerate(optstarlist[varr]) if optmed - num > 2 * optstd]
        for varr2 in optfew:
            if opt_l_threshold[varr][varr2]==minthreshold:
                print(f"few stars in {optcoolist[varr][varr2][:-4]}.fits")

    for varr in infstarlist:
        infmed = statistics.median(infstarlist[varr])
        infstd = statistics.stdev(infstarlist[varr])
        inffew = [i for i, num in enumerate(infstarlist[varr]) if infmed - num > 2 * infstd]
        for varr2 in inffew:
            if inf_l_threshold[varr][varr2]==minthreshold:
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
                #print(f'filename {filename}')
                outfvarr = triangle_match(filename, referencef, outf, match_threshold[varr])
                if outfvarr is None:
                    continue
                opt_matchedf[varr].append(outfvarr)
        else:
            print(f'not execute {varr} trimatch')
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
                if outfvarr is None:
                    continue
                inf_matchedf[varr].append(outfvarr)
        else:
            print(f'not execute {varr} trimatch')
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

    