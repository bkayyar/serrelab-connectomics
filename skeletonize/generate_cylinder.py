"""
Script to generate a binary (b/w) cylinder and store the resulting volume as a numpy file

Takes one command-line argument:
    * Out file: path to the output numpy (.npy) file
"""
import sys
import numpy
from scipy.spatial import distance

CUBE_SIDE = 384

def gen_cylinder(out_file):
    print("Initialzing output volume")
    out_volume = numpy.zeros((CUBE_SIDE, CUBE_SIDE, CUBE_SIDE), dtype='uint8')
    center = numpy.array([CUBE_SIDE/2, CUBE_SIDE/2])
    print("Generating cylinder")
    for z in range(CUBE_SIDE):
        print("Slice " + str(z))
        for y in range(CUBE_SIDE):
            for x in range(CUBE_SIDE):
                pos = numpy.array([x, y])
                if numpy.linalg.norm(center-pos) < 50:
                    out_volume[z, y, x] = 1 
    numpy.save(out_file, out_volume)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise RuntimeError('Usage: python generate_cylinder.py <path_to_output_file>')

    out_file = sys.argv[1]
    gen_cylinder(out_file)

    print("Numpy file written successfully!")
