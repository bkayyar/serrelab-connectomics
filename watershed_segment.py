import sys
import h5py
import numpy
import data_utils

class Thresholds:
    input_path = ""
    output_path = ""
    high_threshold = 0.8
    low_threshold = 0.1
    merge_threshold = 0.2
    merge_size = 800
    dust_size = 600
    is_threshold_relative = True


def baseseg(aff, low_thresh, high_thresh, thresholds, dust_size, is_relative)
    if is_relative:
        low_thresh, high_thresh, thresholds = relative2absolute(aff, low, high, thresholds)
    print("Steepest Ascent")
    seg = steepestascent(aff, low, high)
    print("Divide Plateaus")
    divideplateaus(seg)
    print("Find Basins")
    seg, counts, counts0 = findbasins(seg)
    print("Region Graph")
    rg = regiongraph(aff, seg, len(counts))
    print("Merge Regions")
    new_rg = mergeregions(seg, rg, counts, thresholds, dust_size)
    return seg, new_rg, counts

def atomicseg(thresh):
    h5file = h5py.File(thresh.input_path, 'r')
    raw_data = h5file.get("main")
    aff = data_utils.check_volume(numpy.array(raw_data))
    seg, rg, counts = baseseg(aff, thresh.low_threshold, thresh.high_threshold, 
                                [(thresh.merge_size, thresh.merge_threshold)], 
                                thresh.dust_size, thresh.is_threshold_relative)
    return seg
    


if __name__ == '__main__':
    thresh = Thresholds()
    
    if len(sys.argv) > 2:
        thresh.input_path = sys.argv[1]
        thresh.output_path = sys.argv[2]
        if len(sys.argv) == 8:
            thresh.high_threshold = float(sys.argv[3])
            thresh.low_threshold = float(sys.argv[4])
            thresh.merge_threshold = float(sys.argv[5])
            thresh.merge_size = int(sys.argv[6])
            thresh.dust_size = int(sys.argv[7])
            thresh.is_threshold_relative = False
    else:
        raise RuntimeError('Must pass in at least input and output file names')

    print("==================================")
    print("Input path: " + thresh.input_path)
    print("Output path: " + thresh.output_path)
    print("High Threshold: " + str(thresh.high_threshold))
    print("Low Threshold: " + str(thresh.low_threshold))
    print("Merge Threshold: " + str(thresh.merge_threshold))
    print("Merge Size: " + str(thresh.merge_size))
    print("Dust Size: " + str(thresh.dust_size))
    print("==================================")
        
    
    
