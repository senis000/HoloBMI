

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

from matplotlib import interactive
interactive(True)



def red_channel(red, Afull, new_com, red_im, base_im, fanal, maxdist=3, toplot=True):  
    """
    Function to identify red neurons with components returned by caiman
    red(array-int): mask of red neurons position
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
        plt.savefig(fanal + '/redneurmask.png', bbox_inches="tight")
        
        fig2 = plt.figure()
        R = new_img
        auxA = np.unique(np.arange(Afull.shape[2])*redlabel)
        G = np.transpose(np.nansum(Afull[:,:,auxA],2))
        B = np.zeros((R.shape))
        R = R/np.nanmax(R)
        G = G/np.nanmax(G)
        
        RGB =  np.dstack((R,G,B))
        plt.imshow(RGB)
        plt.savefig(fanal + '/redneurmask_RG.png', bbox_inches="tight")
        plt.close("all") 
            
    return redlabel


def obtain_real_com(Afull, thres=0.1):
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
        new_com[neur,:] = [center_mass[1], center_mass[0]]
                    
    return new_com


def obtain_components(folder, animal, day, estimates, dims):

    # creates folder for plots
    folder_path = os.path.join(folder, animal, day)
    fanal = os.path.join(folder_path, 'plots')
    if not os.path.exists(fanal):
        os.makedirs(fanal)

    # load red.mat
    fmat = os.path.join(folder_path, 'red.mat') 
    redinfo = scipy.io.loadmat(fmat)
    red = redinfo['red']
    base_im = np.reshape(estimates.b, dims).T
    
    # obtain the cnm results
    A_comp = estimates.A[:, estimates.idx_components]
    C_comp = estimates.C_on[estimates.idx_components, :]
    YrA_comp = estimates.YrA[estimates.idx_components, :]

    # convert to normal matrix from sparse
    Afull = np.transpose(np.reshape(A_comp.toarray(),[dims[0], dims[1], A_comp.shape[1]]), [1,0,2])

    # obtain real position of somas (caiman returns center of mass of soma+dendrites)
    new_com = obtain_real_com(Afull)
    
    # import red image
    
    fmat = os.path.join(folder_path, 'red.mat') 
    redinfo = scipy.io.loadmat(fmat)
    red = redinfo['red']
    red_im = redinfo['Im']
    
    # match red neurons with green neurons
    redlabel = red_channel(red, Afull, new_com, red_im, base_im, fanal)  
    ind_red = np.where(redlabel)[0]
    
    # break the sparse matrix to save
    A_comp = scipy.sparse.csr_matrix(np.reshape(Afull.astype('float'), A_comp.shape))
    
    # create dictionary to save as .mat
    f = os.path.join(folder_path, 'redcomp')
    
    dict = {
        'AComp' : A_comp,
        'CComp' : C_comp,   #do not save as int, matlab  does not like it
        'YrA': YrA_comp,
        'com' : new_com.astype('int'),
        'redlabel': redlabel,
        'redIm' : red_im,
        'baseIm' : base_im
        }
    
    redinfo = scipy.io.savemat(f, dict)  

    
def find_index(folder, animal, day, Afull, holofile, com, Ccomp, auxtol=10, cormin=0.6):
    # initialize vars
    neurcor = np.ones((units, all_C.shape[0])) * np.nan
    finalcorr = np.zeros(units)
    finalneur = np.zeros(units)
    finaldist = np.zeros(units)
    iter = 40
	
	# load matlab data
    finfo = os.path.join(folder, animal, day, 'red.mat')  #file name of the mat 
    matinfo = scipy.io.loadmat(finfo)
    redxy = matinfo['red']
    
    fholo = os.path.join(folder, animal, day, holofile)  #file name of the mat 
    holoinfo = scipy.io.loadmat(fholo)
    online_data = holoinfo['holoActivity']
    
	# extract position of neurons
    
    # find distances
    dist = scipy.spatial.distance.cdist(redxy.T, com)
    
    dist[dist>auxtol] = np.nan
    [nmat, ncaim] = np.where(~np.isnan(dist))
    
    for indm, nm in enumerate(np.unique(nmat)):
        possible_match = ncaim[np.where(nmat==nm)[0]]
        neurcor = np.zeros((np.unique(nmat).shape[0], possible_match.shape[0]))
        for indc, nc in enumerate(possible_match):
            auxonline = (np.asarray(online_data[nm, :]) - np.nanmean(online_data[nm, :]))/np.nanmean(online_data[nm, :]) 
            
            neurcor[indm, indc] = pd.DataFrame(np.transpose([auxC[~np.isnan(auxonline)], auxonline[~np.isnan(auxonline)]])).corr()[0][1]
            
	

