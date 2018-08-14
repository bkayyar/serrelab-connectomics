""" Script to convert npy file to GIPL file

Assumptions: Data stored in npy file as float

"""
import sys
import numpy

HEADER_SIZE = 256 #The header is 256 bytes
MAGIC_NUMBER = 4026526128 #A predefined magic number for GIPL files

class GIPLHeader:
    filesize = 0
    sizes = [0.0, 0.0, 0.0, 0.0]
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

def write_metadata(volume, out_file):
    header = GIPLHeader()
    gipl_file = open(out_file, 'wb')
    maxdim = 3

    #Filesize field
    volume_length = numpy.ma.size(volume) #This returns the total number of elements in the volume
    data_type = volume.dtype
    volume_size = data_type.itemsize*volume_length #Number of elements * size of each element
    header.filesize = volume_size + HEADER_SIZE

    #Sizes field
    shape = volume.shape
    assert (len(shape) == 3 or len(shape) == 4)
    if len(shape) == 3:
        maxdim = 3
        for i in range(3):
            header.sizes[i] = shape[i]
        header.sizes[3] = 1.0
    else:
        maxdim = 4
        for i in range(4):
            header.sizes[i] = shape[i]
    for j in range(4):
        gipl_file.write(struct.pack(">H", header.sizes[j]))

    #Image type is float TODO check this 
    header.image_type = 64
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
    gipl_file.write(struct.pack(">f", header.interslice_offset))
    gipl_file.write(struct.pack(">f", header.user_def2))
    gipl_file.write(struct.pack(">I", header.magic_number))
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python npy_to_gipl.py <in_file> <out_file>")
        sys.exit(1)

    in_file = sys.argv[1]
    volume = numpy.load(in_file)
    out_file = sys.argv[2]
    header = write_metadata(volume, out_file)
