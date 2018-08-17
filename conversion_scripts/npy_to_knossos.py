import sys
import numpy
from knossos_utils import knossosdataset as knosdata

if __name__ == '__main__':
    if len(sys.argv) < 4:
        raise RuntimeError('Usage: python npy_to_knossos.py <path_to_input_file> <path_to_output_dataset> <experiment_name>')

    in_file = sys.argv[1]
    out_path = sys.argv[2]
    exp_name = sys.argv[3]
    volume = numpy.load(in_file)
    kd = knosdata.KnossosDataset()
    kd.initialize_from_matrix(out_path, [1.0, 1.0, 1.0], exp_name, data=volume)
