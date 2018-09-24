"""
Script to convert all Numpy npy files in a subdirectory to GIPL files
Command line arguments:
    * --path Path to input npy files 

The wildcard string for glob2 is currently set to glob files from a directory structure 3 deep. We also assume that numpy files end with the .npy extension.
WARNING: GIPL files are said to only support the following data types for our purposes: signed and unsigned 8-bit integers, signed and unsigned 16-bit integers,
signed and unsigned 32-bit integers and 32-bit floats. Make sure your numpy array data type matches one of these!

"""

import os 
# from glob2 import glob
import numpy
import struct
import argparse

HEADER_SIZE = 256 #The header is 256 bytes
MAGIC_NUMBER = 4026526128 #A predefined magic number for GIPL files

class GIPLHeader:
    sizes = [0, 0, 0, 0]
    image_type = 0
    scales = [1.0, 1.0, 1.0, 1.0]
    patient = "No patient information " + " "*57
    matrix = [0.0]*20
    orientation = 0
    par2 = 0
    voxmin = 0.0
    voxmax = 0.0
    origin = [0.0, 0.0, 0.0, 0.0]
    pixval_offset = 0.0
    pixval_cal = 0.0
    interslicegap = 0.0
    user_def2 = 0.0
    magic_number = MAGIC_NUMBER
    volume = []

def find_name(data_type):
    name_type = {"binary": 1, "uint8": 8, "int8": 7, "int16": 15, "uint16": 16} #ITK-SNAP only supports numbers <= 16 bits
    try:
        return name_type[data_type]
    except KeyError:
        raise RuntimeError("Bad numpy data type!")
    
    
def write_file(in_file, out_file, data_type):
    volume = numpy.load(in_file)
    PIXEL_TYPE = find_name(data_type)
    header = GIPLHeader()
    gipl_file = open(out_file, 'wb')
    maxdim = 3

    #Sizes field
    shape = volume.shape
    assert (len(shape) == 3 or len(shape) == 4)
    if len(shape) == 3:
        maxdim = 3
        for i in range(3):
            header.sizes[i] = shape[i]
        header.sizes[3] = 1
    else:
        maxdim = 4
        for i in range(4):
            header.sizes[i] = shape[i]
    for j in range(4):
        gipl_file.write(struct.pack(">H", header.sizes[j]))

    #Image type is stored at the top of the file 
    header.image_type = PIXEL_TYPE
    gipl_file.write(struct.pack(">H", header.image_type))

    #Scales was [2.0, 2.0, 2.0] in Berson's files. TODO find out what this does
    for i in range(maxdim):
        header.scales[i] = 2.0
    for i in range(4):
        gipl_file.write(struct.pack(">f", header.scales[i]))

    #Patient name
    for j in range(80):
        gipl_file.write(struct.pack(">c", header.patient[i]))

    #Matrix. Not sure what this does but this was [0.0]*20
    for i in range(20):
        gipl_file.write(struct.pack(">f", header.matrix[i]))
    
    #Rest of the fields are all 0
    gipl_file.write(struct.pack(">B", header.orientation))
    gipl_file.write(struct.pack(">B", header.par2))
    gipl_file.write(struct.pack(">d", header.voxmin))
    gipl_file.write(struct.pack(">d", header.voxmax))
    for i in range(4):
        gipl_file.write(struct.pack(">d", header.origin[i]))
    gipl_file.write(struct.pack(">f", header.pixval_offset))
    gipl_file.write(struct.pack(">f", header.pixval_cal))
    gipl_file.write(struct.pack(">f", header.interslicegap))
    gipl_file.write(struct.pack(">f", header.user_def2))
    gipl_file.write(struct.pack(">I", header.magic_number))

    gipl_file.seek(HEADER_SIZE, 0) #Seek to end of header
    volume.astype(data_type).transpose(2,1,0).tofile(gipl_file) 
    gipl_file.close()

def main(input_path, wc, data_type):
    print("Globbing files...")
    # files = glob(os.path.join(input_path, wc))
    files = [input_path]
    total_len = len(files)
    count = 0
    for in_file in files:
        out_file = '.'.join(in_file.split('.')[:-1])+".gipl"
        write_file(in_file, out_file, data_type)
        count += 1
        print "Converted  "+in_file+", %d of %d."%(count, total_len)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', dest='input_path',type=str, default=None,
        help='Path to npy segments to be converted')
    parser.add_argument('--wc', dest='wc', type=str, default='**/**/**/*.npy',
        help='Glob2 wildcard to npy segments (** for recursive).')
    parser.add_argument('--dt', dest='data_type',type=str, default=None,
        help='Data type to interpret voxels')
    args = parser.parse_args()
    main(**vars(args))


