from libmetgem.mgf import read as read_mgf
from libmetgem.msp import read as read_msp
from libmetgem.filter import filter_data_multi
from libmetgem.cosine import compute_similarity_matrix

import argparse
import sys
import os
import numpy as np

def get_file_format(input_file):
    ext = os.path.splitext(input_file)[1].lower()
    if ext in ('.mgf', '.msp'):
        return ext[1:]
    
    with open(input_file, 'r') as f:
        head = next(f)
        if head.startswith('BEGIN IONS'):
            return 'mgf'
        elif head.startswith('NAME:'):
            return 'msp'
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MetGem - Molecular Network')
    
    parser.add_argument('-i', '--input', required=True, help='Input file to process')
    parser.add_argument('-o', '--output', default='output.npz', help='Output filename')
    parser.add_argument('-t', '--mz_tolerance', default=0.02, type=float, help='Maximun difference (in Da) between tzo ions masses to consider they correspond to the same ion.')
    
    parser.add_argument('-f', '--use_filtering', default=True, action='store_true', help='Use Filtering')
    parser.add_argument('--min_intensity', type=int, default=0, help='Filter out peaks with relative intensity below this percentage of highest intense peak')
    parser.add_argument('--min_matched_peaks', type=int, default=4, help='Minimum number of common peaks between two spectra')
    parser.add_argument('--parent_filter_tolerance', type=int, default=17, help='in Da')
    parser.add_argument('--min_matched_peaks_search', type=int, default=6, help="Window rank filter's parameters: for each peak in the spectrum, it is kept only if it is in top `min_matched_peaks_search` in the +/-`matched_peaks_window` window")
    parser.add_argument('--matched_peaks_window', type=int, default=50, help='in Da')
    
    parser.add_argument('--is_ms1_data', default=False, action='store_true', help='')
    
    args = parser.parse_args()
   
    mzs = []
    spectra = []

    file_format = get_file_format(args.input)
    
    if file_format == 'mgf':
        read = read_mgf
        mz_keys = ['pepmass']
    elif file_format == 'msp':
        read = read_msp
        mz_keys = ['precursormz', 'exactmass', 'mw']
    else:
        raise ValueError('Unknown file format')
        sys.exit(1)
         
    for params, data in read(args.input, ignore_unknown=True):
        if not args.is_ms1_data:
            mz_parent = 0
            for key in mz_keys:
                try:
                    mz_parent = params[key]
                except KeyError as e:
                    pass
            mzs.append(mz_parent)
        else:
            mzs.append(0)

        spectra.append(data)
        
    if args.use_filtering:
        spectra = filter_data_multi(mzs, spectra, args.min_intensity, args.parent_filter_tolerance,
                                    args.matched_peaks_window, args.min_matched_peaks_search)

    scores = compute_similarity_matrix(mzs, spectra,
                                              args.mz_tolerance, args.min_matched_peaks)
    with open(args.output, 'wb') as f:
        np.savez_compressed(f, scores=scores)
