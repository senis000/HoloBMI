import numpy as np
import matplotlib.pyplot as plt
import cv2
import os

import caiman as cm
from caiman.source_extraction import cnmf as cnmf
from caiman.paths import caiman_datadir


try:
    if __IPYTHON__:
        print("Detected iPython")
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass

from matplotlib import interactive
interactive(True)


def obtain_spatial_filters(folder, fr, use_CNN=True):
    folder_path = folder + 'raw/' + animal + '/' + day + '/'
    fnames = [filename for filename in os.listdir(folder_path) if filename.endswith('.tiff')]  #TODO check .tif or .tiff
    
    # Parameters
    decay_time = .75  # approximate length of transient event in seconds
    gSig = [4, 4]  # expected half size of neurons
    p = 1  # order of AR indicator dynamics
    thresh_CNN_noisy = 0.8  # CNN threshold for candidate components
    gnb = 2  # number of background components
    init_method = 'cnmf'  # initialization method
    
    # set up CNMF initialization parameters

    init_batch = 400  # number of frames for initialization
    patch_size = 96  # size of patch
    stride = 24  # amount of overlap between patches
    K = 25  # max number of components in each patch
    
    SNR_lowest =  2  # very minimum
    min_SNR = 2.5
    rval_lowest = 0.75 #minimum required space correlation.
    rval_thr = 0.8
    cnn_lowest = 0.6 #minimum required CNN threshold.
    min_cnn_thr = 0.9
    
    motion_correct = True # flag for motion correction
    pw_rigid = True
    
    
    
    only_init = True  # whether to run only the initialization

    params_dict = {'fr': fr,
                   'fnames': fname,
                   'decay_time': decay_time,
                   'gSig': gSig,
                   'p': p,
                   'gnv': gnb,
                   'min_SNR': min_SNR,
                   'SNR_lowest': SNR_lowest,
                   'nb': gnb,
                   'init_batch': init_batch,
                   'init_method': init_method,
                   'rf': patch_size//2,
                   'stride': stride,
                   'sniper_mode': True,
                   'thresh_CNN_noisy': thresh_CNN_noisy,
                   'K': K,
                   'rval_lowest': rval_lowest,
                   'rval_thr': rval_thr,
                   'cnn_lowest': cnn_lowest,
                   'min_cnn_thr': min_cnn_thr,
                   'motion_correct':motion_correct
                   }
    opts = cnmf.params.CNMFParams(params_dict=params_dict)
    cnm = cnmf.online_cnmf.OnACID(params=opts)
    cnm.fit_online()
    
    if use_CNN:
        # threshold for CNN classifier
        opts.set('quality', {'min_cnn_thr': min_cnn_thr})
        cnm.estimates.evaluate_components_CNN(opts)
    
    
    