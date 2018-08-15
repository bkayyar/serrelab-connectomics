"""
This program converts image files from the GIPL file format to the HDF5 data format.
Due to a lack of formal documentation of the GIPL file standards, the function that reads GIPL files is a Python transliteration of the Matlab function found a https://www.mathworks.com/matlabcentral/fileexchange/16407-gipl-toolbox
Pass in a GIPL filename as the first argument. The script writes a file called "gipl_converted.hdf5".
"""

import sys
import h5py
import numpy
import struct
from operator import mul

HEADER_SIZE = 256

class GIPLFile:
    filesize = 0
    sizes = []
    image_type = 0
    scales = []
    patient = ""
    matrix = []
    orientation = 0
    par2 = 0
    voxmin = 0.0
    voxmax = 0.0
    origin = []
    pixval_offset = 0.0
    pixval_cal = 0.0
    interslicegap = 0.0
    user_def2 = 0.0
    magic_number = 0

def ReadFileHeader(fname):
    header = GIPLFile()
    OFFSET = 256
    MAGIC_NUMBER = 4026526128

    trans_type = {1:'binary', 7:'char', 8:'uchar', 15:'short', 16:'ushort', 31:'uint',
    32:'int', 64:'float', 65:'double', 144:'C_short', 160:'C_int', 192:'C_float',  
    193:'C_double', 200:'surface', 201:'polygon'}

    trans_orien = {0+1:'UNDEFINED', 1+1:'UNDEFINED_PROJECTION', 2+1:'APP_PROJECTION', 
    3+1:'LATERAL_PROJECTION', 4+1:'OBLIQUE_PROJECTION', 8+1:'UNDEFINED_TOMO', 
    9+1:'AXIAL', 10+1:'CORONAL', 11+1:'SAGITTAL', 12+1:'OBLIQUE_TOMO'}
    
    try:
        gipl_file = open(fname, 'rb')
    except:
        print("Error opening file - check filename maybe?")
        sys.exit(1)

    print("Reading GIPL Header...")

    #Get file size
    gipl_file.seek(0, 2);
    header.filesize = gipl_file.tell()
    gipl_file.seek(0, 0);

    #Sizes header field
    for l in range(4):
        sizes_byte = gipl_file.read(2)
        header.sizes.append(struct.unpack(">H", sizes_byte)[0])
    if (header.sizes[3] == 1):
        maxdim = 3 
    else:
        maxdim = 4
    #header.sizes = header.sizes[:maxdim]

    #Image type header file
    image_type_byte = gipl_file.read(2)
    header.image_type = struct.unpack(">H", image_type_byte)[0]

    #Scales header field
    for m in range(4):
        scales_byte = gipl_file.read(4)
        header.scales.append(struct.unpack(">f", scales_byte)[0])
    header.scales = header.scales[:maxdim]
    
    #Patient header field
    for i in range(80):
        patient_byte = gipl_file.read(1)
        header.patient += struct.unpack(">c", patient_byte)[0]

    #Matrix header field
    for j in range(20):
        matrix_byte = gipl_file.read(4)
        header.matrix.append(struct.unpack(">f", matrix_byte)[0])

    #Orientation header field
    orientation_byte = gipl_file.read(1)
    header.orientation = struct.unpack(">B", orientation_byte)[0]
    
    #par2 header field
    par2_byte = gipl_file.read(1)
    header.par2 = struct.unpack(">B", par2_byte)[0]

    #voxmin and voxmax header fields
    voxmin_byte = gipl_file.read(8)
    header.voxmin = struct.unpack(">d", voxmin_byte)[0]
    voxmax_byte = gipl_file.read(8)
    header.voxmax = struct.unpack(">d", voxmax_byte)[0]

    #Origin header field
    for k in range(4):
        origin_byte = gipl_file.read(8)
        header.origin.append(struct.unpack(">d", origin_byte)[0])
    header.origin = header.origin[:maxdim]

    #Pixval_offset and pixval_cal header fields
    pixval_offset_byte = gipl_file.read(4)
    pixval_cal_byte = gipl_file.read(4)
    header.pixval_offset = struct.unpack(">f", pixval_offset_byte)[0]
    header.pixval_cal = struct.unpack(">f", pixval_cal_byte)[0]

    #Interslice gap header field
    interslice_byte = gipl_file.read(4)
    header.interslice = struct.unpack(">f", interslice_byte)[0]

    #user_def2 header field
    user_def2_byte = gipl_file.read(4)
    header.user_def2 = struct.unpack(">f", user_def2_byte)[0]

    #magic number!!
    magic_number_byte = gipl_file.read(4)
    header.magic_number = struct.unpack(">I", magic_number_byte)[0]
    if (header.magic_number == MAGIC_NUMBER):
        print("Magic Number Matches! Clean Header")
    else:
        raise RunTimeError("MAGIC NUMBER MISMATCH!! BAD FILE")

    gipl_file.close()

    print("=====================================")
    print("filename : " + fname)
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
    print("interslice gap : " + str(header.interslice))
    print("user_def2 : " + str(header.user_def2))
    print("par2 : " + str(header.par2))
    print("Header size (offset): " + str(OFFSET))
    print("=====================================")
    return header

def WriteVolume(header, in_file, out_file):
    try:
        gipl_file = open(in_file, 'rb')
    except:
        print("Error opening file - check filename maybe?")
        sys.exit(1)
    
    format_string = ">"

    if (header.image_type == 1):
        voxelbits = 1
        format_string += "?"
    elif (header.image_type == 7):
        voxelbits = 8
        format_string += "C"
    elif (header.image_type == 8):
        voxelbits = 8
        format_string += "B"
    elif (header.image_type == 15):
        voxelbits = 16
        format_string += "h"
    elif (header.image_type == 16):
        voxelbits = 16
        format_string += "H"
    elif (header.image_type == 31):
        voxelbits = 32
        format_string += "I"
    elif (header.image_type == 32):
        voxelbits = 32
        format_string += "i"
    elif (header.image_type == 64):
        voxelbits = 64
        format_string += "f"
    elif (header.image_type == 65):
        voxelbits = 64
        format_string += "d"

    print("Reading volume...")
    dtype = numpy.dtype(format_string)
    gipl_file.seek(HEADER_SIZE, 0)
    #Try using fromfile or memmap here
    volume = numpy.fromfile(gipl_file, dtype=dtype)
    print("Writing npy file...")
    numpy.save(out_file, volume)
    gipl_file.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python gipl_to_npy.py <path_to_input_file> <path_to_output_file>")
        sys.exit(1)

    in_file = str(sys.argv[1])
    out_file = str(sys.argv[2])
    header = ReadFileHeader(in_file)
    WriteVolume(header, in_file, out_file)
