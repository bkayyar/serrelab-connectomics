"""Script to skeletonize a single 3D object. Takes in two command line arguments:
    * Input filename: This must be a numpy .npy file with a binary (b/w) volume
    * Output filename: This will be another .npy file

Uses scikit's skeletonize_3d function
"""

import sys
import numpy
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize_3d

def skeletonize_volume(in_file, out_file):
    print("Loading volume from input file...")
    input_vol = numpy.load(in_file)
    skeletonized = skeletonize_3d(input_vol)
    numpy.save(out_file, skeletonized)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise RuntimeError("Usage: python skeletonize_segments.py <path_to_input_file> <path_to_output_file>")
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    skeletonize_volume(in_file, out_file)


