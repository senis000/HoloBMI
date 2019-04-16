

import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy
import copy
import tifffile
import cv2
import os

from caiman.components_evaluation import evaluate_components_CNN
from caiman.motion_correction import motion_correct_iteration


interactive(True)



def red_channel(red, nerden, Afull, new_com, red_im, base_im, fanal, maxdist=3, toplot=True):  
    """
    Function to identify red neurons with components returned by caiman
    red(array-int): mask of red neurons position
    nerden(array-bool): array of bool labelling as true components identified as neurons.
    new_com(array): position of the neurons
    red_im(array): matrix MxN: Image of the red channel 
    base_im(array): matrix MxN: Image of the green channel 
    fanal(str): folder where to store the analysis sanity check
    maxdist(int): spatial tolerance to assign red label to a caiman component
    returns
    redlabel(array-bool): boolean vector labelling as True the components that are red neurons 
    """
    #function to identify red neurons
    red_neur = []
    maskred = copy.deepcopy(np.transpose(red)).astype('float32')
    
    # motion correction between the green channel and the red channel 
    # for some reason the motion correction sometimes works one way but not the other
    _, _, shift, _ = motion_correct_iteration(base_im.astype('float32'), red_im.astype('float32'),1)
    
    if np.nansum(abs(np.asarray(shift))) < 20:  # hopefully this motion correction worked
        maskred[:,0] -= shift[1].astype('float32')
        maskred[:,1] -= shift[0].astype('float32')
        # creates a new image with the shifts found
        M = np.float32([[1, 0, -shift[1]], [0, 1, -shift[0]]])
        min_, max_ = np.min(red_im), np.max(red_im)
        new_img = np.clip(cv2.warpAffine(red_im, M, (red_im.shape[0], red_im.shape[1]), flags=cv2.INTER_CUBIC, borderMode=cv2.BORDER_REFLECT), min_, max_)
        print ('Success! shift was: ' + str(shift))
    else:
        print ('Trying other way since shift was: ' + str(shift))
        # do the motion correctio the other way arround
        new_img, _, shift, _ = motion_correct_iteration(red_im.astype('float32'), base_im.astype('float32'),1)
        if np.nansum(abs(np.asarray(shift))) < 20:
            maskred[:,0] += shift[1].astype('float32')
            maskred[:,1] += shift[0].astype('float32')
        else:
            print ('didnt work with shift: ' + str(shift))
            new_img = red_im
            #somehow it didn't work either way
            print ('There was an issue with the motion correction')
        
    
    # find distances
    neur = Afull.shape[2]
    dists = np.zeros((neur,maskred.shape[0]))
    dists = scipy.spatial.distance.cdist(new_com, maskred)
    
    # identfify neurons based on distance
    redlabel = np.zeros(neur).astype('bool')
    aux_redlabel = np.zeros(neur)
    iden_pairs = []  # to debug
    for nn in np.arange(neur):
        if nerden[nn]:
            aux_redlabel[nn] = np.sum(dists[nn,:]<maxdist)
            redn = np.where(dists[nn,:]<maxdist)[0]
            if len(redn):
                iden_pairs.append([nn, redn[0]])  # to debug
    redlabel[aux_redlabel>0]=True
    auxtoplot = new_com[redlabel,:]

    if toplot:
        imgtoplot = np.zeros((new_img.shape[0], new_img.shape[1]))
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1,2,1)
        for ind in np.arange(maskred.shape[0]):
            auxlocx = maskred[ind,1].astype('int')
            auxlocy = maskred[ind,0].astype('int')
            imgtoplot[auxlocx-1:auxlocx+1,auxlocy-1:auxlocy+1] = np.nanmax(new_img)
        ax1.imshow(new_img + imgtoplot, vmax=np.nanmax(new_img))
        
        imgtoplot = np.zeros((new_img.shape[0], new_img.shape[1]))
        ax2 = fig1.add_subplot(1,2,2)
        for ind in np.arange(auxtoplot.shape[0]):
            auxlocx = auxtoplot[ind,1].astype('int')
            auxlocy = auxtoplot[ind,0].astype('int')
            imgtoplot[auxlocx-1:auxlocx+1,auxlocy-1:auxlocy+1] = np.nanmax(new_img)
        ax2.imshow(new_img + imgtoplot, vmax=np.nanmax(new_img))
        plt.savefig(fanal + str(plane) + '/redneurmask.png', bbox_inches="tight")
        
        fig2 = plt.figure()
        R = new_img
        auxA = np.unique(np.arange(Afull.shape[2])*redlabel)
        G = np.transpose(np.nansum(Afull[:,:,auxA],2))
        B = np.zeros((R.shape))
        R = R/np.nanmax(R)
        G = G/np.nanmax(G)
        
        RGB =  np.dstack((R,G,B))
        plt.imshow(RGB)
        plt.savefig(fanal + str(plane) + '/redneurmask_RG.png', bbox_inches="tight")
        plt.close("all") 
            
    return redlabel


def obtain_real_com(Afull, img_size = 20, thres=0.1):
    """
    Function to obtain the "real" position of the neuron regarding the spatial filter
    Afull(array): matrix with all the spatial components
    thres(int): tolerance to identify the soma of the spatial filter
    minsize(int): minimum size of a neuron. Should be change for types of neurons / zoom / spatial resolution
    Returns
    new_com(array): matrix with new position of the neurons
    """
    #function to obtain the real values of com
    new_com = np.zeros((Afull.shape[2], 2))
    for neur in np.arange(Afull.shape[2]):
        center_mass = scipy.ndimage.measurements.center_of_mass(Afull[:,:,neur]>thres)
        if np.nansum(center_mass)==0 :
            center_mass = scipy.ndimage.measurements.center_of_mass(Afull[:,:,neur])
        new_com[neur,:] = [center_mass[0], center_mass[1]]

                    
    return new_com


def obtain_components(folder, A_comp, dims, cnn_tol=0.75):

    # creates folder for plots
    fanal = folder + 'plots/'
    if not os.path.exists(fanal):
        os.makedirs(fanal)

    # load red.mat
    fmat = folder + 'red.mat' 
    redinfo = scipy.io.loadmat(fmat)
    red = redinfo['red']

    # convert to normal matrix from sparse
    Afull = np.reshape(A_comp.toarray(),dims)
    # obtain components from dendrites
    pred, _ = evaluate_components_CNN(A_comp, dims[:2], [4,4])
    nerden = np.zeros(Afull.shape[2]).astype('bool')
    nerden[np.where(pred[:,1]>cnn_tol)] = True
    # obtain real position of somas (caiman returns center of mass of soma+dendrites)
    new_com = obtain_real_com(Afull)
    # match red neurons with green neurons
    redlabel = red_channel(red, nerden, Afull, new_com, red_im, base_im, fanal)  
    ind_red = np.where(redlabel)[0]
    
    # create dictionary to save as .mat
    f = folder + 'redcomp'
    dict = {
        'redLabel' : redlabel,
        'indRed' : ind_red,
        'AComp' : A_comp[:, redlabel],   #do not save as int, matlab  does not like it
        'com' : new_com[redlabel, :].astype('int')
        }
    
    redinfo = scipy.io.savemat(f, dict)  

    
