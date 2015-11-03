#! /usr/bin/python
#
import math
import os
import time
import sys
import numpy as np
from math import *
import pyfits
import pylab as pl
from pylab import *
import matplotlib as plt
############################################################

def chi2(fmodel,fobs,fileout):
    N = len(fobs)
    if N == 1:
        m = 0 # 1 free parameter column density
    else:
        m = 1
        
    tmp,hdr = pyfits.getdata("/home/wwaalkes/radex-models/"+fmodel[0]+".fits", 0,header=True)
    #print "/Users/willwaalkes/Desktop/Radex/data/fits-radex-models-lines/"+fmodel[0]+".fits"
    ncd = tmp.shape[0]
    #print ncd    
    model = np.zeros(shape=(len(fobs),ncd,tmp.shape[1],tmp.shape[2]))
    obs = np.zeros(shape=(len(fobs),tmp.shape[1],tmp.shape[2]))
    #noi = np.zeros(shape=(len(fobs),tmp.shape[1],tmp.shape[2]))
    #snr = np.zeros(shape=(len(fobs),tmp.shape[1],tmp.shape[2]))
    #sig = np.zeros(shape=(len(fobs),tmp.shape[1],tmp.shape[2]))

    for i in range(len(fobs)):
        #print len(fobs)
        #print fobs
        model[i,:,:,:] = pyfits.getdata("/Users/willwaalkes/Desktop/HCN_Research/Radex-Results/"+fmodel[i]+".fits", 0)
        #print model[i,:,:,:]
        #print "/home/wwaalkes/observations/taurus-1-"+fobs[i]+"-area-sig-reproj.fits"
        #if (os.path.isfile("/home/wwaalkes/observations/taurus-1-"+fobs[i]+"-area-sig-reproj.fits")):
        #    obsi = pyfits.getdata("/home/wwaalkes/observations/taurus-1-"+fobs[i]+"-area-sig-reproj.fits", 0)
        #    obs[i,:,:] = obsi.squeeze()
        #    noi[i,:,:] = pyfits.getdata("/home/wwaalkes/observations/taurus-1-"+fobs[i]+"-area-noi-reproj.fits", 0)
        #    snr[i,:,:] = pyfits.getdata("/home/wwaalkes/observations/taurus-1-"+fobs[i]+"-area-snr-reproj.fits", 0)

            #sig[i] = np.mean(obs[0,0:30,0:30])
        #else:
        #    obs[i,:,:] = 0
        #    model[i,:,:,:] = 0
        #    sig[i,:,:] = 1
 
    #sig = np.sqrt(noi**2 + (obs/snr)**2)

    chi2 = np.zeros(shape=(ncd,model.shape[2],model.shape[3]))

    for k in range(ncd): 
        for iobs in range(len(fobs)):
            chi2[k,:,:] += (obs[iobs,:,:] - model[iobs,k,:,:])**2/sig[iobs,:,:]**2
            #print model[iobs,k,:,:]
     #remember to create output folder
    
   # cdmin = 1e10 # SAME AS IN RADEX
   # cdmax = 5e13
   # cd = cdmin*((cdmax/cdmin)**(np.double(np.arange(0,ncd))/ncd))
    cd = np.linspace(1e10,1e13,100)

#    chi2_interp = np.zeros(shape=(500,model.shape[2],model.shape[3]))
        
#==============================================================================
#     #interpolation
#     xp = cd
#     cd_interp = cdmin*((cdmax/cdmin)**(np.double(np.arange(0,500))/500))
#     for i in range(chi2.shape[1]):
#         for j in range(chi2.shape[2]):  
#             yp = chi2[:,i,j]
#             chi2_interp[:,i,j] = np.interp(cd_interp,xp,yp)
#             if((i==30)&(j==50)):   
#             #    print xp
#             #    print cd_interp
#               #  plt.scatter(xp,yp,color='k')  
#               #  plt.scatter(cd_interp,chi2_interp[:,i,j],color='b')
#               #  plt.show()
#             
    chi2_min = np.amin(chi2,axis=0)
    
    chi2_red = chi2_min/(N-m)
    pyfits.writeto("/Users/willwaalkes/Desktop/HCN_Research/N-Best-Fit/"+fileout+"-chi2.fits", chi2_red, hdr, clobber=True)
# 
#     chi2_interp_min = np.amin(chi2_interp,axis=0)
#     pyfits.writeto("/Users/willwaalkes/Desktop/Radex/data/N-best-fit/"+fileout+"-interp-chi2.fits", chi2_interp_min, hdr, clobber=True)
#==============================================================================

    N = np.zeros(shape=(tmp.shape[1],tmp.shape[2]))
#    N_interp = np.zeros(shape=(tmp.shape[1],tmp.shape[2]))
    for i in range(N.shape[0]):
        for j in range(N.shape[1]):
            N[i,j] = cd[np.argmin(chi2[:,i,j])] 
           # N_interp[i,j] = cd_interp[np.argmin(chi2_interp[:,i,j])] 

    pyfits.writeto("/Users/willwaalkes/Desktop/HCN_Research/N-Best-Fit/"+fileout+".fits", N, hdr, clobber=True)
  #  pyfits.writeto("/Users/willwaalkes/Desktop/Radex/data/N-best-fit/"+fileout+"-interp.fits", N_interp, hdr, clobber=True)
 
############################################################

def HCN(sou,source):
    model = [sou+"-radex-hcn-531.7163",sou+"-radex-hcn-620.3040",sou+"-radex-hcn-708.8770",
    sou+"-radex-hcn-797.4333",sou+"-radex-hcn-885.9707",sou+"-radex-hcn-974.4872"]
    obs = [64.068,50.168,33.512,20.475,17.065,8.2064] #in K*km/s
    fileout = source+"-hcn"
    chi2(model,obs,fileout)

#############################################################

def H13CN(sou,source):
    model = [sou+"-radex-h13cn-",sou+"-radex-h13cn-"]
    obs = []
    fileout = source+"-h13cn"
    chi2(model,obs,fileout)
    
#############################################################
    
def HC15N(sou,source):
    model = [sou+"-radex-ch3oh-A-96.7414",sou+"-radex-ch3oh-A-145.1032"]
    obs = []
    fileout = source+"-N-hc15n"
    chi2(model,obs,fileout)
    
#############################################################

#Orion South

HCN("Orion-S","Orion-S")
#H13CN("Orion-S","Orion-S")
#HC15N("Orion-S","Orion-S")

#Orion KL

#HCN("Orion-KL","Orion-KL")
#H13CN("Orion-KL","Orion-KL")
#HC15N("Orion-KL","Orion-KL")

#NGC 6334

#HCN("NGC-6334","NGC-6334")
#H13CN("NGC-6334","NGC-6334")
#HC15N("NGC-6334","NGC-6334")