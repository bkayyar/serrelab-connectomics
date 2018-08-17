import sys
import numpy
from knossos_utils import knossosdataset as knosdata
    
if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise RuntimeError('Usage: python knossos_to_npy.py <path_to_knossos_conf_file> <path_to_output_npy_file>')

    conf_file = sys.argv[1] 
    npy_file = sys.argv[2]

    print("Initializing Knossos dataset...")
    kd = knosdata.KnossosDataset()
    kd.initialize_from_knossos_path(conf_file)
    #Conversion to numpy array probably unnecessary but I'm yet to explore how numpy.save treats non-numpy arrays
    print("Getting volume from raw data...")
    vol = numpy.asarray(kd.from_raw_cubes_to_matrix(size=kd.boundary, offset=[0, 0, 0])) 
    print("Saving numpy file...")
    numpy.save(npy_file, vol)
    print("Numpy file written!")
