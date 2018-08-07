import numpy
import disjoint-sets

""" 
This is a transliteration of the code from files in Seung's Watershed.jl/src/. 
Comments at the top of functions are probably copied straight from the original Julia files
with some necessary modifications and paring down. 
Functions in this file:
steepestascent
divideplateaus
findbasins
"""

"""
The first function attempts to construct the steepest ascent graph from a given affinity graph. 

Inputs:
* `aff`: affinity graph (undirected and weighted). 4D array of affinities, where last dimension is of size 3
* `low`: edges with affinity <= `low` are removed
* `high`: affinities >= `high` are considered infinity

Returns:
* `sag`: steepest ascent graph (directed and unweighted). `sag[x,y,z]` contains 6-bit number enocoding edges outgoing from (x,y,z)

IMPORTANT: Julia is 1-indexed but Python is 0-indexed, so the affinity graph convention probably changes as so.
We follow the convention that:

* `aff[x,y,z,0]` is affinity of voxels at [x-1,y,z] and [x,y,z]
* `aff[x,y,z,1]` is affinity of voxels at [x,y-1,z] and [x,y,z]
* `aff[x,y,z,2]` is affinity of voxels at [x,y,z-1] and [x,y,z]
"""

def steepestascent(aff, low, high):
    assert numpy.ma.size(aff, 4) == 3 #This gives us the size of the matrix along axis 4
    (xdim, ydim, zdim) = aff.shape[:3] #Get the size of the affinity graph (first three axes of aff)
    sag = numpy.zeros((xdim, ydim, zdim), dtype='uint32')

    for z in range(zdim):
        for y in range(ydim):
            for x in range(xdim):
                negx = aff[x, y, z, 0] if x > 0 else low
                negy = aff[x, y, z, 1] if y > 0 else low
                negz = aff[x, y, z, 2] if z > 0 else low
                posx = aff[x+1, y, z, 0] if x < (xdim-1) else low
                posy = aff[x, y+1, z, 1] if y < (ydim-1) else low
                posz = aff[x, y, z+1, 2] if z < (zdim-1) else low
            
                m = numpy.amax([negx, negy, negz, posx, posy, posz])

                if m > low:
                    if (negx == m or negx >= high):
                        sag[x, y, z] |= 0x01
                    if (negy == m or negy >= high):
                        sag[x, y, z] |= 0x02
                    if (negz == m or negz >= high):
                        sag[x, y, z] |= 0x04
                    if (posx == m or posx >= high):
                        sag[x, y, z] |= 0x08
                    if (posy == m or posz >= high):
                        sag[x, y, z] |= 0x10
                    if (posz == m or posz >= high):
                        sag[x, y, z] |= 0x20
    return sag
    
"""
Modify steepest ascent graph so as to
1. Divide non-maximal plateaus into paths that exit as quickly as possible
2. Break ties between multiple outgoing edges

Input: 
`sag`: steepest ascent graph (directed and unweighted). `sag[x,y,z]` contains 6-bit number encoding edges outgoing from (x,y,z)

Returns:
Modified sag
"""
def divideplateaus(sag):
    (xdim, ydim, zdim) = sag.shape()

    dir_array = [-1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim]    
    dir_mask  = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20]
    idir_mask = [0x08, 0x10, 0x20, 0x01, 0x02, 0x04]
    bfs = []

    """ 
    Julia has flat indexing of arrays, which seems to be the Fortran-style "row index changes fastest". I try to replicate this
    by flattening the python array and indexing into it using a scalar, rather than doing some index arithmetic
    """
    flat_sag = numpy.ravel(sag, 'F') #'F' for 'F'ortran
    total_length = len(flat_sag)

    #Queue vertices that have purely outgoing edges
    for idx in range(total_length):
        for d in range(6):
            if (flat_sag[idx] & dir_mask[d]) != 0: #Outgoing edge exists
                if (flat_sag[idx+dir_array[d]] & idir_mask[d]) == 0: #No incoming edge
                    flat_sag[idx] |= 0x40
                    bfs.append(idx)
                    break
    
            
    #Divide plateaus
    bfs_len = len(bfs)
    for bfs_index in range(bfs_len):
        idx = bfs[bfs_index]
        to_set = 0
        for d in range(6):
            if (flat_sag[idx] & dir_mask[d]) != 0: #Outgoing edge exists
                if (flat_sag[idx+dir_array[d]] & idir_mask[d]) != 0: #Incoming edge exists
                    if (flat_sag[idx+dir_array[d] & 0x40) == 0:
                        bfs.append(idx+dir_array[d])
                        flat_sag[idx+dir_array[d]] |= 0x40
                    else:
                        to_set = dir_mask[d]
        flat_sag[idx] = to_set
    sag = numpy.reshape(flat_sag, (xdim, ydim, zdim), 'F')
    return sag

"""
Find basins of attraction

Inputs:
* `sag`: steepest ascent graph (directed and unweighted). `sag[x,y,z]` contains 6-bit number encoding edges outgoing from (x,y,z)

Returns:
* `seg`: segmentation into basins.  Each element of the 3D array contains a *basin ID*, a nonnegative integer ranging from 0 to the number of basins.
* `counts`: number of voxels in each basin
* `counts0`: number of background voxels

A value of 0 in `seg` indicates a background voxel, which has no edges
at all in the steepest ascent graph.  All such singletons are given
the same ID of 0, although they are technically basins by themselves.

The algorithm starts from an unassigned voxel, and identifies all
downstream voxels via breadth-first search (BFS). The search
terminates in two possible ways:

1. downstream voxel that was previously assigned a basin ID =>
assign that ID to queued voxels.
2. no more downstream voxels => assign new basin ID to queued voxels.

Then the queue is emptied, and BFS begins anew at an unassigned voxel.
The algorithm ends when all voxels are assigned.

The 7th bit (0x40) is used to indicate whether a voxel has been
visited during BFS.

The MSB indicates whether a voxel has been assigned a basin ID.  The MSB definition is given in the functions at the end of the file for UInt32 and UInt64 cases.

`findbasins` is applied to the steepest ascent graph after modification by `divideplateaus`  By this point all paths are unique, except in maximal plateaus.
"""
def findbasins(sag):
    (xdim,ydim,zdim) = sag.shape()
    dir_array = [-1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim]
    dir_mask  = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20]
    counts0 = 0  # number of background voxels
    counts = []  # voxel counts for each basin
    bfs = []
    next_id = 1   # initialize basin ID
 
    flat_sag = numpy.ravel(sag, 'F') #To enable flat indexing once again - 'F' for 'F'ortran
    total_length = len(flat_sag)

    for idx in range(total_length):
        if flat_sag[idx] == 0: #No edges - background voxel
            flat_sag[idx] |= 0x80000000 #This is used to get the MSB, under the assumption that sag is composed of uint32. TODO: Check this assumption
            counts0 += 1
        elif (seg[idx] & 0x80000000) == 0:  # not yet assigned
            bfs.append(idx)     # enqueue
            seg[idx] |= 0x40    # mark as visited
            bfs_len = len(bfs)
            for bfs_index in range(bfs_len):
                me = bfs[bfs_index]
                for d in range(6):
                    if (flat_sag[me] & dir_mask[d]) != 0:  # outgoing edge
                        him = me + dir_array[d]  # target of edge
                        if (flat_sag[him] & 0x80000000) != 0: # already assigned
                            for it in bfs:  # assign entire queue to same ID
                                flat_sag[it] = flat_sag[him]  # including high bit
                            counts[flat_sag[him] & 0x80000000] += length(bfs);
                            bfs = []  # empty queue
                            break
                        elif (seg[him] & 0x40) == 0:  # not visited
                            seg[him] |= 0x40    # mark as visited
                            bfs.append(him)    # enqueue
                        # else ignore since visited (try next direction)
                bfs_index += 1      # go to next vertex in queue
            if len(bfs) != 0     # new basin has been created
                counts.append(len(bfs))
                for it in bfs
                    seg[it] = 0x80000000 | next_id    # assign a basin ID
                next_id += 1
                bfs = []

    print("Found: ", str(next_id-1)," components")

    for idx in range(total_length): 
        flat_sag[idx] &= 0x7fffffff     # clear MSB
    end

    seg = numpy.reshape(flat_sag, (xdim, ydim, zdim), 'F')
    return seg, counts, counts0

"""
merge small regions by agglomerative clustering

    new_rg = mergeregions(seg, rg, counts, thresholds, dust_size = 0)

Inputs:
* `seg` - segmentation.  IDs of foreground regions are 1:length(counts).  ID=0 for background.  This is modified in place by the clustering.
* `rg`: region graph as list of edges, array of (weight,id1,id2) tuples. The edges should be presorted so that weights are in descending order. Region IDs should be consistent with those in `seg`, except no zeros.
* `counts`: sizes of regions in `seg` (modified in place)
* `thresholds`: sequence of (size_th,weight_th) pairs to be used for merging
* `dust_size`: after merging, tiny regions less than dust_size to be eliminated by changing them to background voxels

Returns:
* `new_rg`: new region graph after clustering, same format as `rg`.

Agglomerative clustering proceeds by considering the edges of the region graph in sequence.  If either region has size less than `size_th`, then merge the regions. When the weight of the edge in the region graph is less than or equal to `weight_th`, agglomeration proceeds to the next `(size_th,weight_th)` in `thresholds` or terminates if there is none.
"""
def mergeregions(seg, rg, counts, thresholds, dust_size=0):
    counts_len = len(counts)
    sets = DisjointSets(counts_len) #Using DisjointSets object from third-party code. Must check TODO
    for (size_th, weight_th) in thresholds:
        for (weight, id1, id2) in rg:
            s1 = sets.find(id1)
            s2 = sets.find(id2)
            if (weight > weight_th) and (s1 != s2):
                if (counts[s1] < size_th) or (counts[s2] < size_th):
                    counts[s1] += counts[s2]
                    counts[s2] = 0
                    sets.union(s1, s2)
                    s = sets.find(s1)   # this is either s1 or s2
                    (counts[s], counts[s1]) = (counts[s1], counts[s]) #Is this just swapping elements?
    print("Done merging")

    # define mapping from parents to new segment IDs
    # and apply to redefine counts
    remaps = numpy.zeros(counts_len)     
    next_id = 1 
    for idx in range(counts_len):
        s = sets.find(idx)
        if (remaps[s] == 0) and (counts[s] >= dust_size): 
            remaps[s] = next_id
            counts[next_id] = counts[s]    
            next_id += 1 
    counts = counts[:next_id-1] #TODO ensure this is the same as resize!()

    # apply remapping to voxels in seg
    # note that dust regions will get assigned to background
    flat_seg = numpy.ravel(seg, 'F')
    for idx in range(len(flat_seg)):
        if flat_seg[idx] != 0:    # only foreground voxels
            flat_seg[idx] = remaps[sets.find(flat_seg[idx])]
    seg = numpy.reshape(flat_seg, seg.shape(), 'F')
    print("Done with remapping, total: ", str(next_id-1), " regions")

    # apply remapping to region graph
    in_rg = [[0]] * (next_id-1)
    new_rg = []
    for (weight, id1, id2) in rg:
        s1 = remaps[sets.find(id1)]
        s2 = remaps[sets.find(id2)]
        if (s1 != s2) and (s1 != 0) and (s2 != 0):  # ignore dust regions
            (s1,s2) = (min(s1,s2), max(s1, s2))
            if s2 not in in_rg[s1]:
                new_rg.append((weight, s1, s2))
                in_rg[s1].append(s2)
    print("Done with updating the region graph, size: ", str(len(new_rg)))
    return new_rg



