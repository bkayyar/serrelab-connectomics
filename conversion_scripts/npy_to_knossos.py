"""
Script to convert Numpy npy files to Knossos datasets (.conf and .raw files).
Command line arguments:
    * Path to input numpy file
    * Path to output dataset. Knossos creates top-level directories (mag1, mag2, etc.) automatically - provide the parent of this top-level directory! Do not provide uncreated directory paths.
    * Experiment Name

WARNING: I don't think using npy_to_knossos.py and knossos_to_npy.py in conjunction will leave you with identical npy files - 
once touched by Knossos, forever changed (Knossos uses uint8 I think). That is to say, npy1 -> knossos (.raw) -> npy2 =/=> npy1 = npy2
if you use npy_to_knossos.py and knossos_to_npy.py
"""
import sys
import numpy
from knossos_utils import knossosdataset as knosdata

if __name__ == '__main__':
    if len(sys.argv) < 4:
        raise RuntimeError('Usage: python npy_to_knossos.py <path_to_input_file> <path_to_output_dataset> <experiment_name>')

    in_file = sys.argv[1]
    out_path = sys.argv[2]
    exp_name = sys.argv[3]
    volume = numpy.load(in_file).T
    kd = knosdata.KnossosDataset()
    kd.initialize_from_matrix(out_path, [1.0, 1.0, 1.0], exp_name, data=volume)
