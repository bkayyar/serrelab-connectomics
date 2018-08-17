# serrelab-connectomics

Example commands to use conversion scripts. Look at the comments in each file for argument descriptions.

Knossos to Numpy:
python knossos_to_npy.py ./datasets/mag1/knossos.conf ./datasets/npy_files/exp_top.npy

Numpy to Knossos:
python npy_to_knossos.py ./datasets/npy_files/exp_numpy.npy ./datasets/ First-Experiment

NOTE: The npy_to_knossos.py script automatically creates the top-level Knossos dataset directory. _*Do Not*_ supply the
name of an uncreated top-level directory.
