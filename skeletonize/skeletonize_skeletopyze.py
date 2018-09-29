"""Script to skeletonize segments. Takes in two command line arguments:
    * Input filename: This must be a numpy .npy file
    * Output filename: This will be another .npy file

Uses the skeletopyze library 
"""

from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy
import matplotlib.pyplot as plt
import skeletopyze

MEMBRANE_ID = 2081

def skeletonize_volume(in_file, out_file):

    fig = plt.figure()
    ax = fig.gca(projection="3d")
    fig2 = plt.figure()
    ax2 = fig2.gca(projection="3d")
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    print("Loading segments from input file...")
    input_segments = numpy.load(in_file).astype("int8")
    ax2.voxels(input_segments, edgecolor='k')
    params = skeletopyze.Parameters()
    skel_coords = skeletopyze.get_skeleton_graph(input_segments, params)
    for n in skel_coords.nodes():
        ax.scatter(skel_coords.locations(n).x(), skel_coords.locations(n).y(), skel_coords.locations(n).z())
        #ax.plot(skel_coords.locations(n).x(), skel_coords.locations(n).y(), skel_coords.locations(n).z())
        #print("(%d, %d, %d)"%(skel_coords.locations(n).x(), skel_coords.locations(n).y(), skel_coords.locations(n).z()))

    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise RuntimeError("Usage: python skeletonize_segments.py <path_to_input_file> <path_to_output_file>")
    
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    skeletonize_volume(in_file, out_file)


