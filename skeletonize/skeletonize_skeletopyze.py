"""Script to skeletonize segments. Takes in two command line arguments:
    * Input filename: This must be a numpy .npy file
    * Output filename: This will be another .npy file

Uses the skeletopyze library 
"""

import sys
import numpy
import matplotlib.pyplot as plt
import skeletopyze

MEMBRANE_ID = 2081

def skeletonize_volume(in_file, out_file):
    print("Loading segments from input file...")
    input_segments = numpy.load(in_file)
    (zdim, ydim, xdim) = input_segments.shape
    single_segment = numpy.zeros((zdim, ydim, xdim), dtype='uint8') #Arrays to store results
    print("Getting segment IDs...")
    seg_ids = numpy.unique(input_segments) #Get segment IDs
    num_segs = len(seg_ids)
    params = skeletopyze.Parameters()
    for idx in seg_ids:
        if idx != MEMBRANE_ID: #Skip cell membrane
            print("Skeletonizing segment ID " + str(idx) + " out of " + str(num_segs) + "...")
            chosen_seg = idx == input_segments #Get elements for the current segment
            for z in range(zdim):
                for y in range(ydim):
                    for x in range(xdim):
                        if chosen_seg[z, y, x]:
                            single_segment[z, y, x] = 1 #Convert to binary image
            skel_coords = skeletopyze.get_skeleton_graph(single_segment, params)
            #plt.imshow(single_segment[12, :, :], cmap='gray')
            #plt.show()
    #plt.imshow(skeletonized[12, :, :], cmap='gray')
    #plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise RuntimeError("Usage: python skeletonize_segments.py <path_to_input_file> <path_to_output_file>")
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    skeletonize_volume(in_file, out_file)


