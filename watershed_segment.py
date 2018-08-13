import sys
import h5py
import numpy
import data_utils
from graph_functions import *

DEFAULT_LOW = 0.1
DEFAULT_HIGH = 0.8
DEFAULT_DUST_SIZE = 600
DEFAULT_MERGE_SIZE = 800
DEFAULT_MERGE_THRESHOLD = 0.2
DEFAULT_THRESHOLD_RELATIVE = True

class Thresholds:
    input_path = ""
    high_threshold = DEFAULT_HIGH 
    low_threshold = DEFAULT_LOW
    merge_threshold = DEFAULT_MERGE_THRESHOLD
    merge_size = DEFAULT_MERGE_SIZE
    dust_size = DEFAULT_DUST_SIZE
    is_threshold_relative = DEFAULT_THRESHOLD_RELATIVE

def relative2absolute(aff, low, high, thresholds):
    if numpy.ma.size(aff) > 3*1024*1024*128:
        hist = numpy.histogram(aff[1:min(1024, numpy.ma.size(aff, 1)),
                                   1:min(1024, numpy.ma.size(aff, 2)),
                                   1:min(128,  numpy.ma.size(aff, 3)),
                                    :], nbins=10000000)
    else:
        hist = numpy.histogram(aff, nbins=1000000)
    low = percent2thd(h, low)
    high = percent2thd(h, high)
    for i in range(len(thresholds)):
        thresholds[i] = (thresholds[i][0], percent2thd(h, thresholds[i][1]))
    return low, high, thresholds
    

def baseseg(aff, low_thresh, high_thresh, thresholds, dust_size, is_relative)
    if is_relative:
        low_thresh, high_thresh, thresholds = relative2absolute(aff, low, high, thresholds)
    print("Steepest Ascent")
    seg = steepestascent(aff, low, high)
    print("Divide Plateaus")
    seg = divideplateaus(seg)
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
    aff = numpy.array(raw_data)
    h5file.close()
    seg, rg, counts = baseseg(aff, thresh.low_threshold, thresh.high_threshold, 
                                [(thresh.merge_size, thresh.merge_threshold)], 
                                thresh.dust_size, thresh.is_threshold_relative)
    return seg

def watershed(aff, low=DEFAULT_LOW, high=DEFAULT_HIGH, thresholds=[(DEFAULT_MERGE_SIZE, DEFAULT_MERGE_THRESHOLD)], 
              dust_size=DEFAULT_DUST_SIZE, is_threshold_relative=DEFAULT_IS_THRESHOLD_RELATIVE):
    seg, rg, counts = baseseg(aff, low, high, thresholds, dust_size, is_threshold_relative)
    return seg


def mergerg(seg, rg, thd=0.5):
    pd = {}
    num = 0
    for (a,c,p) in rg: #affinity, child, parent
        assert (c>0 and p>0)
        if a >= thd:
            num += 1
        pd[c] = (p,a)
    print("Total number of merging edges: " + str(num))
    print("Total number: " + str(numpy.ma.size(rg)))
    
    rd = {}
    rset = set()
    for (a, c0, p0) in rg:  
        c = c0
        p = p0
        while (a>=thd and p in pd):
            a = pd[p][1]
            p = pd[p][0]
        if p != p0:
            rd[c0] = p
            rset.add(P)
    print("Number of trees: " + str(len(rset)))

    num = 0
    flat_seg = numpy.ravel(seg, 'F')
    flat_seg_len = len(flat_seg)
    # set the segment id as relative root id
    for i in range(flat_seg_len):
        # segment id
        sid = flat_seg[i]
        if sid in rd:
            flat_seg[i] = rd[sid]
            num = num + 1
    print("Really merged edges: " + str(num))

def wsseg2d(affs, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rt=0.5)
    seg = numpy.zeros(size(affs)[:3], dtype='uint32')
    for z in range(numpy.ma.size(affs,2)):
        seg[:,:,z], rg = watershed(affs[:,:,z,:], low, high, thresholds, dust_size)
        seg[:,:,z] = mergerg(seg[:,:,z], rg, thd_rg)
    return seg

def wsseg(affs, dim=3, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rg=0.5)
    assert (dim==2 or dim==3)
    if dim==2:
        return wsseg2d(affs, low, high, thresholds, dust_size, thd_rg)
    else:
        seg, rg = watershed(affs, low, high, thresholds, dust_size)
        seg = mergerg(seg, rg, thd_rg)
        return seg

def rg2segmentpairs(rg):
    flat_rg = numpy.ravel(rg, 'F')
    flat_rg_len = len(flat_rg)
    segmentPairAffinities = numpy.zeros(flat_rg_len, dtype='float32')
    segmentPairs = numpy.zeros((flat_rg_len, 2), dtype='uint32')

    for i in range(flat_rg_len):
        t = rg[i]
        segmentPairAffinities = t[0]
        segmentPairs[i, 0] = t[1]
        segmentPairs[i, 1] = t[2]

    return segmentPairs, segmentPairAffinities
    

if __name__ == '__main__':
    thresh = Thresholds()
    
    if len(sys.argv) > 2:
        thresh.input_path = sys.argv[1]
        output_path = sys.argv[2]
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
    print("Output path: " + output_path)
    print("High Threshold: " + str(thresh.high_threshold))
    print("Low Threshold: " + str(thresh.low_threshold))
    print("Merge Threshold: " + str(thresh.merge_threshold))
    print("Merge Size: " + str(thresh.merge_size))
    print("Dust Size: " + str(thresh.dust_size))
    print("==================================")

    seg = atomicseg(thresh)
    out_file = h5py.File(output_path, 'w')
    out_file.create_dataset('main', data=seg)
    out_file.close()
    print("Output written to file")
