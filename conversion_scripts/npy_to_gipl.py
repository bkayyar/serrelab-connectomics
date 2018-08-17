"""
Script to convert Numpy npy files to GIPL files
Command line arguments:
    * Path to input npy file 
    * Path to output GIPL file. 

WARNING: GIPL files are said to only support the following data types for our purposes: signed and unsigned 8-bit integers, signed and unsigned 16-bit integers,
signed and unsigned 32-bit integers and 32-bit floats. Make sure your numpy array data type matches one of these!
"""

import sys
import numpy
import struct

PIXEL_TYPE = 0 #Data types are determined later
DTYPE_STRING = '' 
HEADER_SIZE = 256 #The header is 256 bytes
MAGIC_NUMBER = 4026526128 #A predefined magic number for GIPL files

class GIPLHeader:
    filesize = 0
    sizes = [0, 0, 0, 0]
    image_type = 0
    scales = [1.0, 1.0, 1.0, 1.0]
    patient = ""
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

def find_name(dtype):
    name_type = {"uint8": 8, "int8": 7, "int16": 15, "uint16": 16, "uint32": 31, "int32": 32, "float32": 64}
    name = dtype.name
    try:
        return name_type[name]
    except KeyError:
        raise RuntimeError("Bad numpy data type!")
    
    
def write_file(volume, out_file):
    header = GIPLHeader()
    print("Opening file...")
    gipl_file = open(out_file, 'wb')
    maxdim = 3

    print("Calculating header fields...")
    #Filesize field
    volume_length = numpy.ma.size(volume) #This returns the total number of elements in the volume
    data_type = numpy.dtype(DTYPE_STRING)
    volume_size = data_type.itemsize*volume_length #Number of elements * size of each element
    header.filesize = volume_size + HEADER_SIZE

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
    header.patient = "No patient information " + " "*57
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

    print("Writing volume...")
    gipl_file.seek(HEADER_SIZE, 0) #Seek to end of header
    #'unsafe' reminds you that all types of data conversion are allowed. Transpose to resolve Fortran-vs-C array semantics
    volume.astype(DTYPE_STRING, casting='unsafe').T.tofile(gipl_file) 

    gipl_file.close()
    print("File written successfully!")
    return header
    
if __name__ == '__main__':

    trans_type = {1:'binary', 7:'char', 8:'uchar', 15:'short', 16:'ushort', 31:'uint',
    32:'int', 64:'float', 65:'double', 144:'C_short', 160:'C_int', 192:'C_float',
    193:'C_double', 200:'surface', 201:'polygon'}

    trans_orien = {0+1:'UNDEFINED', 1+1:'UNDEFINED_PROJECTION', 2+1:'APP_PROJECTION',
    3+1:'LATERAL_PROJECTION', 4+1:'OBLIQUE_PROJECTION', 8+1:'UNDEFINED_TOMO',
    9+1:'AXIAL', 10+1:'CORONAL', 11+1:'SAGITTAL', 12+1:'OBLIQUE_TOMO'}

    if len(sys.argv) < 3:
        print("Usage: python npy_to_gipl.py <in_file> <out_file>")
        sys.exit(1)

    in_file = sys.argv[1]
    volume = numpy.load(in_file)
    DTYPE_STRING = str(volume.dtype)
    PIXEL_TYPE = find_name(volume.dtype) 
    out_file = sys.argv[2]
    header = write_file(volume, out_file)

    print("=====================================")
    print("output filename : " + out_file)
    print("filesize : " + str(header.filesize))
    print("sizes : " + str(header.sizes)[:])
    print("image_type : " + str(header.image_type) + '-' + trans_type[header.image_type])
    print("Scales : " + str(header.scales)[:])
    print("patient : " + header.patient)
    print("matrix : " + str(header.matrix)[:])
    print("orientation : " + str(header.orientation) + '-' + trans_orien[header.orientation+1])
    print("voxel_min : " + str(header.voxmin))
    print("voxel_max : " + str(header.voxmax))
    print("origin : " + str(header.origin)[:])
    print("pixval_offset : " + str(header.pixval_offset))
    print("pixval_cal : " + str(header.pixval_cal))
    print("interslice gap : " + str(header.interslicegap))
    print("user_def2 : " + str(header.user_def2))
    print("par2 : " + str(header.par2))
    print("Header size (offset): " + str(HEADER_SIZE))
    print("=====================================")

