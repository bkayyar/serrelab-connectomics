"""Script to skeletonize segments. Takes in two command line arguments:
    * Input filename: This must be a numpy .npy file
    * Output filename: This will be another .npy file

Uses scikit's skeletonize_3d function

We assume that the segment ID that occurs most frequently in the volume belongs to 
the segment that has the largest volume.
"""

import sys
import numpy
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize_3d

NUM_SEGS = 5 #Number of largest segments to skeletonize

def skeletonize_volume(in_file, out_file):
    print("Loading segments from input file...")
    input_segments = numpy.load(in_file)
    (zdim, ydim, xdim) = input_segments.shape
    single_segment = numpy.zeros((zdim, ydim, xdim), dtype='uint8') #Arrays to store results
    skeletonized = numpy.zeros((zdim, ydim, xdim), dtype='uint8')
    segments = numpy.zeros((zdim, ydim, xdim), dtype='uint8')
    print("Getting most frequent segment IDs...")
    (seg_ids, counts) = numpy.unique(input_segments, return_counts=True) #Get segment IDs and their counts
    sorted_count_indices = numpy.argpartition(counts, -NUM_SEGS)[-NUM_SEGS:] #This returns the seg_ids array indices of the IDs of the NUM_SEGS largest segments 
    largest_segs = seg_ids[sorted_count_indices] #Get the seg IDs of the NUM_SEGS largest segments
    descending_segs = largest_segs[numpy.argsort(-largest_segs)] #Sort in descending order to get the largest segment first
    count = 0
    for idx in descending_segs:
        if idx != 2098: #Skip 0 class - cell membrane
            count += 1
            print("Skeletonizing segment ID " + str(idx) + ", " + str(count) + " out of " + str(NUM_SEGS) + "...")
            for z in range(zdim):
                for y in range(ydim):
                    for x in range(xdim):
                        if input_segments[z, y, x] == idx:
                            single_segment[z, y, x] = 1 #Convert to binary image
            single_segment_skeletonized = skeletonize_3d(single_segment) #Skeletonize single segment
            segments = numpy.add(segments, single_segment) #Merge entire segments, not skeletonized
            skeletonized = numpy.add(single_segment_skeletonized, skeletonized) #Merge current segment skeleton with other segment skeletons
            single_segment = numpy.zeros((zdim, ydim, xdim), dtype='uint8') #Reinitialize single segment array
            #plt.imshow(single_segment[12, :, :], cmap='gray')
            #plt.show()
    numpy.save(out_file, skeletonized)
    numpy.save("full_segments.npy", segments) #Save full segments without skeletonization
    #plt.imshow(skeletonized[12, :, :], cmap='gray')
    #plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise RuntimeError("Usage: python skeletonize_segments.py <path_to_input_file> <path_to_output_file>")
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    skeletonize_volume(in_file, out_file)


